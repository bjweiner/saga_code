
; takes a list of objects with Hecto spectrum filenames,
; and fitting each line in turn.
; return an array of arrays of fit params and errors

pro do_linefit_loop, objid, z, prefx, specfile, aperture, wrest, fitpars, fiterrors, fitquals

  nspec = n_elements(specfile)
  wfit = wrest
  sighecto = 2.0
  cdummy = ' '
  q = 0
;  fitpars = []
;  fiterrors = []
  fitquals = intarr(nspec)

  for i = 0,nspec-1 do begin
     jj = aperture[i]-1
     fname = prefx + specfile[i] 
     lam = mrdfits(fname,0,/silent) 
     flux = mrdfits(fname,1,/silent) 
     ivar = mrdfits(fname,2,/silent) 
     mask = mrdfits(fname,4,/silent) * 0.0 
     pixels = where(mask[*,jj] eq 0) 
     lam1 = lam[pixels,jj] 
     lam0 = lam1 / (1.0 + z[i]) 
     flux1 = flux[pixels,jj] 
     err1 = 1.0/sqrt(ivar[pixels,jj]) 
     fitline_generic, lam1, flux1, wfit, fitpars1, fiterrors1,  $
                      fluxerr=err1, zfit=z[i], siginst=sighecto, $
                      chisq=chisq, dof=dof, /plotlabels, /multiline, $
                      /freenii
     print, objid[i]," EWs = ",fitpars1[4]/fitpars1[0],fitpars1[5]/fitpars1[0] 
     if i eq 0 then begin
        fitpars = fitpars1
        fiterrors = fiterrors1
     endif else begin
        fitpars = [fitpars, fitpars1]
        fiterrors = [fiterrors, fiterrors1]
     endelse
     read, q, prompt="Enter quality grade: "
     fitquals[i] = q
;     read, cdummy, prompt="Hit enter for next: "
     
  endfor

end


  
