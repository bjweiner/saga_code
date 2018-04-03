
; fit a line in a 1-d spectrum with known redshift

; parameters are 
; continuum, slope, location, sigma, intensity

; fitline_generic is a variant of fitline1d.pro that is not grism
; specific. Rather than passing in a structure, just pass in 
; arrays of wavelength and flux (and a redshift in the parameter zfit)
; an instrumental resolution siginst in A, and optionally an error
; vector in fluxerr

pro fitline_generic, lambda, flux, wrest, fitpars, fiterrors, fluxerr = fluxerr, zfit = zfit, siginst = siginst, fixwave = fixwave, fixsigma = fixsigma, plotlabels = plotlabels, fluxcorr = fluxcorr, scalederrors = scalederrors, chisq = chisq, dof = dof, multiline = multiline, trimdata = trimdata, nspec1d = nspec1d, freenii = freenii

  forward_function gauss_slopecont, gauss_oiii, gauss_ha

  defplotcolors
  if not keyword_set(plotlabels) then plotlabels = 0
  if not keyword_set(fluxcorr) then fluxcorr = 1.0
; Can scale data, errors, and initial guesses up by some
; factor.  This may help give more reliable errors from mpfit ?
;   fluxscale = 1.0e18
  fluxscale = 1.0

  if not keyword_set(fluxerr) then begin
     fluxerr = flux-flux + 1.0
  endif

  if keyword_set(multiline) and multiline ne 0 then begin
     if wrest gt 5000. and wrest lt 5010. then begin
        linefunc = 'gauss_oiii'
        ifunc = 1
        npars = 6
        fitradlo = 200.0
        fitradhi = 100.0
     endif else if wrest gt 6559. and wrest lt 6566. then begin
        linefunc = 'gauss_ha'
; gauss_ha_nii allows NII to be free in fit, SII is optional
        if keyword_set(freenii) then linefunc = 'gauss_ha_nii'
        ifunc = 2
        npars = 6
        fitradlo = 100.0
        fitradhi = 200.0
     endif else begin
        linefunc = 'gauss_slopecont'
        ifunc = 0
        npars = 5
        fitradlo = 200.0
        fitradhi = 200.0
     endelse
  endif else begin
     linefunc = 'gauss_slopecont'
     ifunc = 0
     npars = 5
     fitradlo = 200.0
     fitradhi = 200.0
  endelse        
  if keyword_set(trimdata) then begin
     fitradlo = fitradlo * trimdata
     fitradhi = fitradhi * trimdata
  endif

  wlabel = [ 3728.0, 3933., 3967., 4340., 4861.3, 5006.8, 6562.8 ]
  nlabel = n_elements(wlabel)
  wlabelstr = string(fix(wlabel))

;  if not keyword_set(zfit) then zfit = gtarg.z
  if not keyword_set(zfit) then zfit = 0.0
  wobs = wrest * (1+zfit)
; moderately conservative limits of good data in spectrum
;  wmin = 10900.0
;  wmax = 16600.0
; allow more or less any wavelength of spectrum
  wmin = 100.0
  wmax = 1.0e5
; sigma of line in A due to instr resolution and object size
  if not keyword_set(siginst) then siginst = 1.4
  gsig = siginst
; radius of window in which to use data to fit the line
; depends on if we are fitting more than one line
;  fitrad = 0.08e4

; which 1-d spectrum to use
;  if not keyword_set(nspec1d) then nspec1d = 0
;  if nspec1d+1 gt gtarg.nspectra1d then nspec1d = gtarg.nspectra1d-1
;  if nspec1d+1 gt gtarg.nspectra1d then nspec1d = 0
  
; how much data to plot.  trim the region used to determine
; y plot limits since spectrum often turns up at ends.
;  iplot = where(gtarg.spec1dlambda[nspec1d,*] gt wmin and gtarg.spec1dlambda[nspec1d,*] lt wmax, ntoplot)
;  iplotlim = where(gtarg.spec1dlambda[nspec1d,*] gt wmin+300. and gtarg.spec1dlambda[nspec1d,*] lt wmax-300., ntoplotlim)
;  iplot = where(lambda gt wmin and lambda lt wmax, ntoplot)
;  iplotlim = where(lambda gt wmin and lambda lt wmax, ntoplotlim)

  wfitmin = max([wmin, wobs-fitradlo])
  wfitmax = min([wmax, wobs+fitradhi])

;  iplot = where(lambda gt wfitmin and lambda lt wfitmax, ntoplot)
  iplot = where(lambda gt wmin and lambda lt wmax, ntoplot)
;  iplotlim = where(lambda gt wfitmin-100 and lambda lt wfitmax+100, ntoplotlim)
  iplotlim = iplot

;  indx = where(gtarg.spec1dlambda[nspec1d,*] gt wfitmin and gtarg.spec1dlambda[nspec1d,*] lt wfitmax, ntofit)
  indx = where(lambda gt wfitmin and lambda lt wfitmax, ntofit)

  if ntofit lt 5 then begin
;     print, "Warning: not enough points for fit at wl ",wobs," object ",gtarg.objname
     print, "Warning: not enough points for fit at wl ",wobs
     fitpars = replicate(0.0, npars)
     fiterrors = replicate(1000.0, npars)
     scalederrors = replicate(1000.0, npars)
     chisq = 1000.0
     dof = 1.0
     nofit = 1
; need to make the plot even if not fitting so keep going
;     return
  endif else begin
     nofit = 0
     wdata = lambda
     fdata = flux 
     edata = fluxerr

; correct for the fluxes being too small by avg exp time ~ 1328 sec?
;  fdata = fdata * 1327.9
;  edata = edata * 1327.9


; use mpfit to fit an emission line
; scale up the fluxes to correct for axe's normalization problem
; if needed (in later reductions, generally unnecessary)
; or initial guesses appropriate for uncorrected input fluxes
     if keyword_set(fluxcorr) then begin
        fdata = fdata * fluxcorr
        edata = edata * fluxcorr
; This was wrong - since I'm scaling the data, no need to scale
; the initial guesses.
;        start = [ 1.0e-21*fluxcorr, 0.0, wobs, gsig, 5.0e-20*fluxcorr ]
     endif
; scale data to make fitting function not have to work on tiny numbers
     fdata = fdata * fluxscale
     edata = edata * fluxscale
;     contguess = 1.0e-18
;     intguess = 5.0e-17
;     contguess = 20.0
;     intguess = 10.0
     contguess = total(fdata[indx]) / n_elements(fdata[indx]) / fluxscale
     intguess = contguess * gsig * 2.0
     if npars eq 5 then begin
        start = [ contguess*fluxscale, 0.0, wobs, gsig, intguess*fluxscale ]
     endif else if npars eq 6 then begin
        start = [ contguess*fluxscale, 0.0, wobs, gsig, intguess*fluxscale, 0.15 ]
     endif else begin
; So far npars should be 5 or 6 and this shouldn't happen
        start = [ contguess*fluxscale, 0.0, wobs, gsig, intguess*fluxscale ]
     endelse

; need this structure to constrain parameters
     parfix = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]}, npars )
; Only allow positive line intensities
     parfix[4].limited = [1,0]
     parfix[4].limits  = [0.D,1.0d10]
     if npars eq 6 then begin
        parfix[5].limited = [1,0]
        parfix[5].limits  = [0.D,1.0d10]
     endif
; to fix the wavelength to input value
; could also consider constraining it to some region
     if keyword_set(fixwave) then parfix[2].fixed = 1
; to fix the sigma
     if keyword_set(fixsigma) then parfix[3].fixed = 1

; return fitted params in fitpars - 
; continuum, slope, location, sigma, intensity
     fitpars = mpfitfun(linefunc, wdata[indx], fdata[indx], edata[indx], start, $
                        parinfo=parfix, dof=dof, bestnorm=fbestnorm, $
                        perror=fiterrors, /quiet)
; if the fit fails (like no data) then it may fail to return
; fiterrors
     if n_elements(fiterrors) le 0 then begin
        fitpars = replicate(0.0, npars)
        fiterrors = replicate(1000.0, npars)
        scalederrors = replicate(1000.0, npars)
        chisq = 1000.0
        dof = 1.0
        nofit = 1
     endif
; If you think the reduced chisq should be ~1 you can scale
; the reported errors to get the "real" errors
     chisq = fbestnorm
     scalederrors = fiterrors * sqrt(chisq / dof)

; end the else-part where fitting is done
  endelse

; make a plot
;  ploterr, wdata, fdata, edata
;  ploterr, gtarg.spec1dlambda[iplot], gtarg.spec1dflux[iplot], $
;        gtarg.spec1derr[iplot]
  toplabel = 'object ' + " z= " + string(zfit)
; Auto-scaling including the error bar messes up if one point has
; a large error
;  ymin = min(gtarg.spec1dflux[iplot]-gtarg.spec1derr[iplot]) * fluxcorr
;  ymax = max(gtarg.spec1dflux[iplot]+gtarg.spec1derr[iplot]) * fluxcorr
  if ntoplot gt 0 then begin
     marg = max(flux[iplotlim]) - min(flux[iplotlim])
     ymin = min(flux[iplotlim]-0.05*marg) * fluxcorr * fluxscale
     ymax = max(flux[iplotlim]+0.05*marg) * fluxcorr *fluxscale
     yplotrange = [ymin-0.05*(ymax-ymin), ymax+0.1*(ymax-ymin)] 
     plot, lambda[iplot], flux[iplot] * fluxcorr*fluxscale, $
           xtitle="wavelength", ytitle="flux, spec1d ", title=toplabel, psym=4, $
           yrange = yplotrange
     errplot, lambda[iplot], $
              (flux[iplot] - fluxerr[iplot]) * fluxcorr*fluxscale, $
              (flux[iplot] + fluxerr[iplot]) * fluxcorr*fluxscale
; if we made a fit, plot the fitted function
     if nofit eq 0 then oplot, wdata, call_function(linefunc, wdata, fitpars), color=!red, thick=2.0
     if plotlabels ne 0 then begin
        ylabel = ymax-0.06*(ymax-ymin)
        for i = 0, nlabel - 1 do begin
           xyouts, wlabel[i]*(1.0+zfit), ylabel, wlabelstr[i], charsize=0.98, $
                   orientation = 90
        endfor
     endif
  endif else begin
; Make an empty plot if there wasn't any data
     plot, [0.0, 1.0], [0.0, 1.0], $
           xtitle="wavelength", ytitle="flux", title=toplabel, psym=4
     xyouts, 0.5, 0.5, "no 1-d data for object ", $
             charsize=2.1
  endelse

     
; correct for flux problem?
;  fitpars[4] = fitpars[4] * 1327.9

end


; function passed to mpfitfun
; parameters of the gaussian are location (central wl), sigma, intensity
; so params of mygauss are continuum, location, sigma, intensity

function mygauss, x, p

  return, p[0] + GAUSS1(x, p[1:3])

end

; Gaussian with linear sloping continuum
; parameters of the gaussian are location, sigma, intensity
; params of gauss_slopecont are
; continuum, cont-slope, location, sigma, intensity
; linear continuum is zero-pointed at the central location

function gauss_slopecont, x, p

  return, p[0] + p[1]*(x-p[2]) + GAUSS1(x, p[2:4])

end





