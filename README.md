# Solar Data Mining

Python tools for mining and cross-matching multi-wavelength solar flare data for event dates of interest. The repository supports retrieval and inspection of radio, EUV, and high-energy solar datasets used in solar flare and space weather studies.

## Purpose

This repository contains analysis utilities developed for compiling event-based solar flare datasets. The tools are designed to support studies connecting solar flare properties with radio bursts, energetic particles, and CME-related space weather outcomes.

## Main capabilities

- Query and process event-date-based multi-wavelength solar observations.
- Retrieve and inspect decameter–hectometric radio dynamic spectra.
- Download SDO/AIA image cutouts for flare events.
- Download and plot NoRP GHz light curves.
- Match MWA observations with HEK flare/event information.
- Generate MWA observing catalogs.
- Support RHESSI flare catalog based hard X-ray event searches.

## Repository contents

| File | Purpose |
|---|---|
| `Multiwb_datamine.py` | Main multi-wavelength solar data mining utilities. |
| `get_radioDS_fordates.py` | Tools for retrieving radio dynamic spectra for selected dates. |
| `downAIAcutouts.py` | Downloads SDO/AIA cutouts for solar events. |
| `downNorP_GHzLC_plot.py` | Downloads and plots NoRP GHz light curves. |
| `MWA_HEK_match.py` | Matches MWA observations with HEK solar event information. |
| `MWA_obscat_maker.py` | Creates MWA observing catalog files. |
| `hessi_flare_catalog.csv` | RHESSI flare catalog table used by some utilities. |
| `time_conv.py` | Helper routines for time conversion. |

## Scientific context

These tools were developed for event-based solar flare and space weather studies, especially for constructing catalogs that connect flare emission, radio bursts, CMEs, and energetic particle events. They are useful for Sun-as-a-star analyses where multi-wavelength flare parameters must be extracted consistently across many events.

## Notes

The scripts are research utilities rather than a packaged software library. Users should inspect file paths, date formats, and data-server dependencies before running them.
