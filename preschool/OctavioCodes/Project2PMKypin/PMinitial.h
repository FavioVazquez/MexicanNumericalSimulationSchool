      PARAMETER (Lblock=4)
      PARAMETER (Nblocks=NROW/Lblock)
      PARAMETER (NSPEC = NROW/2)   ! No waves shorter than Ny
      PARAMETER (marr  = NROW+1)   ! Arrays for FFT
      PARAMETER (mf67  = NROW/2)
      COMMON / MASK/	iMask(Nblocks,Nblocks,Nblocks) ! Refinement Mask
      COMMON / GRID/     GRX(NROW,NROW,NROW),   ! Spectrum of Perurbations
     +                   GRY(NROW,NROW,NROW),
     +                   GRZ(NROW,NROW,NROW)
      COMMON / BUFF /	XPt(Nmaxpart),YPt(Nmaxpart),ZPt(Nmaxpart),
     +			VXt(Nmaxpart),VYt(Nmaxpart),VZt(Nmaxpart)
      COMMON / LevelsM/ Max_lev_abs,Nmin_lev,Nmax_lev







