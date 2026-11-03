'''
This code helps to regularise any matrix given XxY matrix, X array, Y array, (xstart, xend)
& (ystart,yend)
It is assumed that the in put is a matrix in which case one dimension either Row or Column
will keep its size throughout data. 
'''

import numpy as np

def Reg(DS,R,C,r0,r1,c0,c1,reg_ax='both',stats='min'):
	'''
	Reg(DS,R,C,x0,x1,y0,y1)
	
	INPUTS
	------
	
	DS: 2D array of data
	R: numpy array of Row ticks
	C: numpy array of Column ticks
	r0: start Row value required
	r1: end Row value required
	c0: Start Column value required
	c1: End Column value required 	
	reg_ax: Which axes to be regularised. Default is 'both'.
			'both': Both row and column will be uniformly regridded
			'col': only column axis will be regridded assuming Row is fine
			'row': Only row will be regridded assuming column is fine.	
	stats: 'median','min','mean' Which metric should be used on the Row or column array to select the stepsize
			
	Output
	------
	DSn: New regularised DS
	R: New Row. Bad row locations will be masked. Old if its column to be regularised
	C: New column. Bad column locations will be masked. Old if its Row to be regularised
		
	'''
	DCs=C[1:]-C[:-1]
	DRs=R[1:]-R[:-1]

	if stats=='min':
		Stat=np.min
	if stats=='median':
		Stat=np.median
	if stats=='mean':
		Stat=np.mean
				
	if reg_ax=='row':
		gdr=Stat(DRs)
		GRs=np.arange(r0,r1+gdr,gdr)
		if GRs[-1]>r1:
			GRs=GRs[:-1]
		DSn=np.zeros((len(GRs),DS.shape[1]))+np.nan
		gdr=gdr.astype(float)
		cnt=0
		Rmask=[]
		for i in GRs:		
			val=np.min(np.abs((R-i).astype(float)))
			if val>gdr:
				print('Data gap found!! Filled with NaNs')
				cnt+=1
				Rmask+=[True]
				continue
			loc=np.argmin(np.abs((R-i).astype(float)))
			DSn[cnt,:]=DS[loc,:]
			Rmask+=[False]
			cnt+=1
		Rgs=np.ma.array(GRs,mask=Rmask)
		Cgs=np.ma.array(C,mask=False)
		
	if reg_ax=='col':
		gdc=Stat(DCs)
		GCs=np.arange(c0,c1+gdc,gdc)
		if GCs[-1]>c1:
			GCs=GCs[:-1]
		DSn=np.zeros((DS.shape[0],len(GCs)))+np.nan
		gdc=gdc.astype(float)
		cnt=0
		Cmask=[]
		for i in GCs:		
			val=np.min(np.abs((C-i).astype(float)))
			if val>gdc:
				print('Data gap found!! Filled with NaNs')
				cnt+=1
				Cmask+=[True]
				continue
			loc=np.argmin(np.abs((C-i).astype(float)))
			DSn[:,cnt]=DS[:,loc]
			Cmask+=[False]
			cnt+=1
		Cgs=np.ma.array(GCs,mask=Cmask)
		Rgs=np.ma.array(R,mask=False)
		
	if reg_ax=='both':
		gdc=Stat(DCs)
		GCs=np.arange(c0,c1+gdc,gdc)
		if GCs[-1]>c1:
			GCs=GCs[:-1]

		gdr=Stat(DRs)
		GRs=np.arange(r0,r1+gdr,gdr)
		if GRs[-1]>r1:
			GRs=GRs[:-1]

		DSn=np.zeros((len(GRs),len(GCs)))+np.nan
		gdc=gdc.astype(float)
		gdr=gdr.astype(float)
		
		Cmask=[False]*len(GCs)
		Rmask=[False]*len(GRs)
		cntC=0
		for i in GCs:		
			valC=np.min(np.abs((C-i).astype(float)))
			if valC>gdc:
				print('Data gap found!! Filled with NaNs')
				Cmask[cntC]=True
				cntC+=1
				continue	
			locC=np.argmin(np.abs((C-i).astype(float)))
			cntR=0
			for j in GRs:				
				valR=np.min(np.abs((R-j).astype(float)))
				if valR>gdr:
					print('Data gap found!! Filled with NaNs')
					Rmask[cntR]=True
					cntR+=1
					continue
				locR=np.argmin(np.abs((R-j).astype(float)))
				DSn[cntR,cntC]=DS[locR,locC]
				cntR+=1
			cntC+=1	
		Cgs=np.ma.array(GCs,mask=Cmask)
		Rgs=np.ma.array(GRs,mask=Rmask)	
			
	return DSn,Rgs,Cgs