#!/usr/bin/env python
"""
Solar Orbiter / Metis quick-look movie maker.

Downloads Metis L2 visible-light tB files from the Solar Orbiter Archive (SOAR) for
sT..eT, keeps one frame per dT, aligns them (north up, Sun centred), saves a FITS cube,
writes PNG frames and stitches them into an MP4.

How to use: edit the INPUT section below, then run   python metis_quickmovie.py
Needs: numpy, scipy, astropy, matplotlib, ffmpeg (or `pip install imageio-ffmpeg`)
"""
import glob, os, shutil, subprocess, sys, warnings
import urllib.parse, urllib.request
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.time import Time
import astropy.units as u
from scipy.ndimage import map_coordinates

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", RuntimeWarning)
TAP = "https://soar.esac.esa.int/soar-sl-tap/"
RSUN_M = 6.957e8

# =========================================================================================
#                                       INPUT
# =========================================================================================
# --- time range and cadence ---
sT       = "2025-05-30T07:00:00"   # start time (UT)
eT       = "2025-05-31T00:00:00"   # end time (UT)
dT       = 30                      # cadence [min]; frame nearest each step is used

# --- movie type ---
mode     = "base_diff"             # "intensity" | "running_diff" | "base_diff"
nbase    = 3                       # base_diff: no. of first frames median-combined as base
rscale   = 2.0                     # diff images multiplied by r**rscale (0 = off)
clip     = 1.0                     # colour-scale percentile clip [%]

# --- output ---
datadir  = "/Metis/metis_tb"            # where FITS files are downloaded
outdir   = "Metis"          # cube, PNG frames and mp4 go here
fps      = 5                       # movie frames per second
dpi      = 150                     # PNG resolution
# =========================================================================================
# ================================================================ 1. find + download
def query_soar(t0, t1, product="metis-vl-tb"):
    q = (f"SELECT data_item_id, filename, begin_time FROM v_sc_data_item "
         f"WHERE instrument='METIS' AND descriptor='{product}' AND level='L2' "
         f"AND begin_time>='{t0}' AND begin_time<='{t1}' ORDER BY begin_time")
    url = TAP + "tap/sync?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "CSV", "QUERY": q})
    lines = urllib.request.urlopen(url, timeout=120).read().decode().strip().split("\n")[1:]
    return [dict(zip(("item", "file", "time"), l.split(","))) for l in lines if l]


def pick_cadence(rows, t0, t1, dt_min):
    """One file per dT step: the file nearest each grid time, within dT/2."""
    if not rows:
        return []
    tt = Time([r["time"] for r in rows]).mjd
    grid = np.arange(Time(t0).mjd, Time(t1).mjd + 1e-9, dt_min / 1440)
    keep = []
    for g in grid:
        i = int(np.argmin(np.abs(tt - g)))
        if np.abs(tt[i] - g) <= dt_min / 2880 and i not in keep:
            keep.append(i)
    return [rows[i] for i in keep]


def download(rows, outdir):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for k, r in enumerate(rows, 1):
        out = os.path.join(outdir, r["file"])
        if not os.path.exists(out):
            url = TAP + "data?" + urllib.parse.urlencode(
                {"retrieval_type": "PRODUCT", "data_item_id": r["item"], "product_type": "SCIENCE"})
            urllib.request.urlretrieve(url, out + ".part")
            os.replace(out + ".part", out)
        print(f"  [{k}/{len(rows)}] {r['file']}")
        paths.append(out)
    return paths


# ================================================================ 2. read + align
def read_frame(path):
    with fits.open(path) as hdul:
        h = next(x for x in hdul if x.data is not None and x.data.ndim == 2)
        data, hdr = h.data.astype(np.float32), h.header
    data[~np.isfinite(data) | (data <= 0)] = np.nan       # occulter / invalid
    t = Time(next(hdr[k] for k in ("DATE-AVG", "DATE-OBS", "DATE-BEG") if k in hdr))
    rsun = float(hdr["RSUN_ARC"]) if "RSUN_ARC" in hdr else \
        np.degrees(np.arcsin(RSUN_M / float(hdr["DSUN_OBS"]))) * 3600
    return dict(data=data, t=t, rsun=rsun, wcs=WCS(hdr).celestial)


def to_grid(fr, axis):
    """Resample onto a north-up grid in R_sun centred on the Sun (removes roll/drift)."""
    k = u.Unit(fr["wcs"].wcs.cunit[0]).to(u.arcsec)
    X, Y = np.meshgrid(axis, axis)
    px, py = fr["wcs"].wcs_world2pix((X * fr["rsun"] / k).ravel(), (Y * fr["rsun"] / k).ravel(), 0)
    c = np.vstack([py, px])
    img = map_coordinates(np.nan_to_num(fr["data"]), c, order=1).reshape(X.shape)
    ok = map_coordinates(np.isfinite(fr["data"]).astype(float), c, order=0).reshape(X.shape)
    img[ok < 0.5] = np.nan
    return img


def build_cube(paths):
    frames = sorted((read_frame(p) for p in paths), key=lambda f: f["t"].mjd)
    f0 = frames[0]; ny, nx = f0["data"].shape
    k = u.Unit(f0["wcs"].wcs.cunit[0]).to(u.arcsec)
    lon, lat = f0["wcs"].wcs_pix2world([0, nx - 1, 0, nx - 1], [0, 0, ny - 1, ny - 1], 0)
    lon = (lon * k + 648000) % 1296000 - 648000           # wrap to +-180 deg (in arcsec)
    half = np.max(np.abs(np.r_[lon, lat * k])) / f0["rsun"]
    axis = np.linspace(-half, half, int(np.maximum(nx, ny)))
    cube = np.stack([to_grid(f, axis) for f in frames]).astype(np.float32)
    return cube, axis, Time([f["t"] for f in frames])


def save_cube(path, cube, axis, times):
    h = fits.Header()
    for i in (1, 2):
        h[f"CTYPE{i}"], h[f"CUNIT{i}"] = ("SOLAR_X", "SOLAR_Y")[i - 1], "R_SUN"
        h[f"CRPIX{i}"], h[f"CRVAL{i}"], h[f"CDELT{i}"] = 1, axis[0], axis[1] - axis[0]
    tab = fits.BinTableHDU.from_columns([fits.Column("time", "A23", array=times.isot)], name="TIMES")
    fits.HDUList([fits.PrimaryHDU(cube, h), tab]).writeto(path, overwrite=True)


# ================================================================ 3. movie
def make_movie(cube, axis, times, mode, out, fps=5, nbase=3, rscale=2.0, clip=1.0, dpi=150):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = np.hypot(*np.meshgrid(axis, axis))
    if mode == "intensity":
        disp, idx = np.log10(cube), np.arange(len(cube))
    elif mode == "base_diff":
        disp, idx = (cube - np.nanmedian(cube[:nbase], axis=0)) * r ** rscale, np.arange(len(cube))
    elif mode == "running_diff":
        disp, idx = (cube[1:] - cube[:-1]) * r ** rscale, np.arange(1, len(cube))
    else:
        sys.exit("mode must be intensity, running_diff or base_diff")

    v = disp[np.isfinite(disp)]
    if mode == "intensity":                                  # LASCO-C2-style colour table
        vmin, vmax = np.percentile(v, [clip, 100 - clip])
        try:
            import sunpy.visualization.colormaps  # noqa
            cmap = matplotlib.colormaps["soholasco2"]
        except Exception:
            cmap = matplotlib.colormaps["gist_heat"]
        cblabel = r"log$_{10}$ tB"
    else:                                                    # differences: grey, symmetric
        lim = np.percentile(np.abs(v), 100 - clip); vmin, vmax = -lim, lim
        cmap = matplotlib.colormaps["gray"]
        cblabel = "difference" + (rf" $\times r^{{{rscale:g}}}$" if rscale else "")
    cmap = cmap.copy(); cmap.set_bad("black")

    framedir = os.path.splitext(out)[0] + "_frames"
    shutil.rmtree(framedir, ignore_errors=True); os.makedirs(framedir)
    fig = plt.figure(figsize=(7.8, 6.6))
    ax = fig.add_axes([0.15, 0.12, 0.66, 0.84]); cax = fig.add_axes([0.83, 0.12, 0.03, 0.84])
    im = ax.imshow(disp[0], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[axis[0], axis[-1], axis[0], axis[-1]], interpolation="bilinear")
    ax.add_patch(plt.Circle((0, 0), 1, fill=False, color="white", lw=1))
    ax.set_xlabel(r"Solar X ($R_\odot$)", fontsize=17); ax.set_ylabel(r"Solar Y ($R_\odot$)", fontsize=17)
    ax.tick_params(labelsize=15)
    cb = fig.colorbar(im, cax=cax); cb.set_label(cblabel, fontsize=15); cb.ax.tick_params(labelsize=13)
    name = {"intensity": "tB", "base_diff": "base diff", "running_diff": "running diff"}[mode]
    txt = ax.text(0.03, 0.97, "", transform=ax.transAxes, va="top", color="white", fontsize=14,
                  bbox=dict(facecolor="black", alpha=0.6, edgecolor="none"))
    for n, i in enumerate(idx):
        im.set_data(disp[n])
        txt.set_text(f"Metis VL {name}  {times[i].strftime('%Y-%m-%d %H:%M')} UT")
        fig.savefig(os.path.join(framedir, f"frame_{n:04d}.png"), dpi=dpi)
    plt.close(fig)

    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        try:
            import imageio_ffmpeg; ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
        except ImportError:
            sys.exit(f"PNG frames are in {framedir}/ but ffmpeg was not found to make the mp4")
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", os.path.join(framedir, "frame_%04d.png"),
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264",
                    "-preset", "veryfast", "-crf", "18", "-pix_fmt", "yuv420p", out], check=True)
    print(f"movie: {out}  ({len(idx)} frames, PNGs in {framedir}/)")


# ================================================================ main
def main():
    t0, t1 = Time(sT).isot, Time(eT).isot
    print(f"searching SOAR: Metis L2 tB  {t0} -> {t1}")
    rows = pick_cadence(query_soar(t0, t1), t0, t1, dT)
    if len(rows) < 2:
        sys.exit(f"only {len(rows)} frame(s) found for this period/cadence")
    print(f"{len(rows)} frames at ~{dT:g} min cadence, downloading to {datadir}/")
    paths = download(rows, datadir)

    os.makedirs(outdir, exist_ok=True)
    cube, axis, times = build_cube(paths)
    tag = os.path.join(outdir, f"metis_{Time(t0).strftime('%Y%m%dT%H%M')}_{Time(t1).strftime('%Y%m%dT%H%M')}")
    save_cube(tag + "_cube.fits", cube, axis, times)
    print(f"cube: {tag}_cube.fits  {cube.shape}")
    make_movie(cube, axis, times, mode, f"{tag}_{mode}.mp4",
               fps=fps, nbase=nbase, rscale=rscale, clip=clip, dpi=dpi)


if __name__ == "__main__":
    main()
