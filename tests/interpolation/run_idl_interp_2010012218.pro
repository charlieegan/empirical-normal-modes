;
; Operates in four modes:
; 
; 1) Reads analysed theta, winds and PV on eta levels.
; Interpolates vertically to theta levels.
; Calculates isentropic vorticity by finite difference of winds.
; Calculates Montgomery potential by integration of temperature.
; Calculates isentropic density from theta-derivatives of M.
; Calculates PV from vorticity/density and compares with eta coord estimate.
; Finds area, mass and circulation enclosed by PV contours on
; isentropic surfaces.
; Outputs results in bs* files to be read by PV inverter to
; determine the modified Lagrangian mean state by iteration.
;
; Note: PV is modified using function of theta LAIT2PV (in SI units).
;
; 2) Reads analysed theta, winds and PV on eta levels and calculates
; sigma and PV from Montgomery potential as in mode 1.
; Reads integrals from bs* files output in mode 1.
; Reads equivalent latitudes output by PV inverter.
; Reads zonally symmetric basic state output by PV inverter. 
; Calculates wave activity diagnostics and outputs file. 
; Can switch background between day 0 and time evolving (lbsevol).
;
; 3) As for mode 2, but uses output from PV inverter for neutral
; reference state.
;
; 4) Performs both modes 1 and 2 and spawns the PV inverter script 
; inbetween. Avoids re-calculation of the isentropic-coord fields.
;
; Modified from backint.pro
;
; John Methven   v1.1 12/1/05.
; John Methven   v2.1 5/2/07.  (adding mode 2)
; John Methven   v2.2 17/8/07. (adding full p-mom and p-energy calcs)
; John Methven   v2.3 3/5/10.  (adding WA flux calcs)
; John Methven   v3.1 14/9/10. (adding mode 4)
; Hannah Croad        12/02/25 (adding new OT background state)
;
; *************************************************************************
;
; Read tracer header output by tradv program.
;
pro read_trhead, hunit, fhead, itratyp, alon, alat, aeta, beta

openr,hunit,fhead

astr=' '
readf,hunit,astr
print,astr
readf,hunit,astr
readf,hunit,astr,ntrac,format='(A22,I4)'
print,astr,ntrac
itratyp=lonarr(ntrac)
readf,hunit,astr,itratyp,format='(A18,99I4)'
; taurel=fltarr(ntrac)
; readf,hunit,astr,taurel,format='(A30,10F8.3)'
; print,astr,taurel
readf,hunit,astr ; fix for old and new taurel headers
readf,hunit,mg
alon=dindgen(mg)*360./float(mg)
readf,hunit,jgg
alat=dblarr(jgg)
readf,hunit,alat,format='(8E13.5)'
readf,hunit,nlp
aeta=dblarr(nlp)
beta=dblarr(nlp)
readf,hunit,aeta,format='(8E13.5)'
readf,hunit,beta,format='(8E13.5)'

return
end
; *************************************************************************
pro read_trdata,idatadate,fpath,fstemh,fstemo,ltracer $
               ,itratyp,longitude,latitude,aeta,beta,psurf,data

fhead=fstemh+strtrim(idatadate,2)
finput=fstemo+strtrim(idatadate,2)
lr4out=1

hunit=1
scrap=' '
read_trhead,hunit,fpath+fhead,itratyp,longitude,latitude,aeta,beta
ntrac=n_elements(itratyp)
nlon=n_elements(longitude)
nlat=n_elements(latitude)
nlev=n_elements(aeta)-1
;
; Read row compression data from tracer header only.
;
iflag=lonarr(nlat,nlev)
jflag=lonarr(nlat,nlev,ntrac)
if ltracer eq 1 then begin
   readf,hunit,scrap
   for n=1,ntrac do begin
      readf,hunit,iflag,format='(70I1)'
      jflag(*,*,n-1)=iflag(*,*)
   endfor
endif
close,hunit
;
; If surface pressure is diagnosed, save it separately.
;
ntracin=ntrac
itratypin=itratyp
lpsurf=where(itratyp eq 0)
if lpsurf(0) ne -1 then begin
   print,'Diagnostic record ',lpsurf(0),' is surface pressure.'
   itratyp=itratypin(where(itratypin ne 0))
   ntrac=n_elements(itratyp)
endif

print, 'reading ERA data: ', fpath+finput
openr, unit, fpath+finput, /get_lun

iflag=0
if lr4out eq 1 then begin
   fdata=fltarr(nlon)
   data=fltarr(nlon,nlat,nlev,ntrac)
   psurf=fltarr(nlon,nlat)
endif else begin
   fdata=dblarr(nlon)
   data=dblarr(nlon,nlat,nlev,ntrac)
   psurf=dblarr(nlon,nlat)
endelse
psurf(*,*)=1.
data(*,*,*,*)=0.

;
; Note that the swap_endian lines are not required if the TRADV
; output was created on the same computer.
; It is needed if the computer producing the diout* files had a 
; processor using different convention for fp number binary storage.
;
ntr=0
if ltracer eq 1 then begin
   for n=0,ntracin-1 do begin
      for l=0, nlev-1 do begin
	 for j=0,nlat-1 do begin
            if jflag(j,l,n) eq 1 then begin
               readu, unit, fdata
;               fdata=swap_endian(fdata,/swap_if_little_endian)
               data(*,j,l,ntr) = fdata
            endif
         endfor	
      endfor
      ntr=ntr+1
   endfor
endif else begin
   for n=0,ntracin-1 do begin
      case itratypin(n) of
      0: begin
         for j=0,nlat-1 do begin
            readu, unit, fdata
;            fdata=swap_endian(fdata,/swap_if_little_endian)
            psurf(*,j) = fdata
         endfor	
      end
      else: begin
         for l=0, nlev-1 do begin
	    for j=0,nlat-1 do begin
               readu, unit, fdata
;               fdata=swap_endian(fdata,/swap_if_little_endian)
               data(*,j,l,ntr) = fdata
            endfor	
         endfor
         ntr=ntr+1
      end
      endcase
   endfor
endelse
free_lun, unit
;
; If data is global force it to be NH only.
;
lnhonly=0
if lnhonly eq 1 and latitude(nlat-1) lt 0. then begin
    nlatin=nlat
    nlat=nlatin/2
    latitude=latitude(0:nlat-1)
    psurf=psurf(*,0:nlat-1)
    data=data(*,0:nlat-1,*,*)
endif

return
end

; *************************************************************************
pro plotxy,idatadate,data,longitude,latitude,alev $
   ,lsurf,lisen,thlev $
   ,attrchoice,units,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
;
; Procedure to plot horizontal sections through 3D field on maps.
;
nlon=n_elements(longitude)
nlat=n_elements(latitude)
dlon=longitude(1)-longitude(0)
xarr=fltarr(nlon+1)
xarr(0:nlon-1)=longitude(*)
xarr(nlon)=xarr(0)
xydata=fltarr(nlon+1,nlat)
minval=min(data)

if lsurf eq 1 then begin
    xydata(0:nlon-1,*)=data(*,*)
endif else begin
    xydata(0:nlon-1,*)=data(*,*,etachoice)
endelse
xydata(nlon,*)=xydata(0,*)

if llog eq 1 then begin
   minval=-998.
   xycop=xydata
   xycop(*,*)=minval-1.
   nonzero=where(xydata gt 0.)
   if nonzero(0) gt -1 then xycop(nonzero)=alog10(xydata(nonzero))
   xydata=xycop
endif

DDEG=15.
NLONG=360./DDEG
NLATG=180./DDEG
ALONG=FINDGEN(NLONG+1)*DDEG
ALATG=FINDGEN(NLATG)*DDEG-89.999
;
; Note: clevs here is offset by -dlev relative to array
; used for sections due to IDL difference between colours used
; for shading map projections and normal graphs.
; The grid area is surrounded by a border, one point wide,
; with values chosen to close contours properly. 
;
clevarr=lev0+(findgen(nclev)-1)*dlev
linearr=2*(clevarr lt 0.)
zel=where(clevarr eq 0.)
zeroc=zel(0)
clevs=clevarr
lines=linearr
cthick=linearr
cthick(*)=1
if zeroc ne -1 then begin
   if nozer eq 1 then begin
      clevs=fltarr(nclev-1)
      lines=fltarr(nclev-1)
      if zeroc gt 0 then begin
         clevs(0:zeroc-1)=clevarr(0:zeroc-1)
	 lines(0:zeroc-1)=clevarr(0:zeroc-1)
      endif
      if zeroc lt nclev-1 then begin	
	 clevs(zeroc:*)=clevarr(zeroc+1:*)
	 lines(zeroc:*)=linearr(zeroc+1:*)
      endif
   endif else begin
      clevs=clevarr
      lines=linearr
      lines(zeroc)=1
   endelse
endif
if llog eq 1 then lines(*)=0

labs=clevs
labs(*)=0
scal=240./(max(clevs)-min(clevs))

if lisen eq 1 then begin
   thc=thlev(etachoice)
   if thc lt 999. then begin
      thf=strmid(strtrim(thlev(etachoice),2),0,3)
   endif else begin
      thf=strmid(strtrim(thlev(etachoice),2),0,4)
   endelse
   maintit=units+'   !7h!x = ' $
       +thf+'K   '+strtrim(idatadate,2)
endif else begin
   maintit=units+'   !7g!x = ' $
       +strtrim(alev(etachoice),2)+'   '+strtrim(idatadate,2)
endelse

case projchoice of
0: begin
   map_set, 0, 0, /cylindrical, title=maintit
end
1: begin
   limplot=[35, -180, 90, 180]
   map_set, 90, 0, /stereographic, /isotropic $
          , title=maintit, limit=limplot
end
2: begin
   map_set, -90, 0, /stereographic, /isotropic, title=maintit
end
3: begin
   seglon=dlon*nlon
   lowlat=15.
   xarr1seg=xarr
   xydata1seg=xydata
   xarr=fltarr(2*nlon+1)
   xydata=fltarr(2*nlon+1,nlat)
   xarr(0:nlon-1)=xarr1seg(0:nlon-1)-seglon
   xarr(nlon:2*nlon-1)=xarr1seg(0:nlon-1)
   xarr(2*nlon)=seglon
   xydata(0:nlon-1,*)=xydata1seg(0:nlon-1,*)
   xydata(nlon:2*nlon,*)=xydata1seg(*,*)
   nlong=2*seglon/ddeg
   along=-seglon+ddeg*findgen(nlong+1)
   limplot=[lowlat, -seglon, 90, seglon]
   map_set, 90, 0, /stereographic, /isotropic, /noerase $
          , title=maintit $
          , /noborder, limit=limplot
   th1=where(clevs eq 280.)
   if th1(0) ne -1 then begin
      cthick(th1(0))=3
   endif
end
endcase

if lfill gt 0 then begin
   tmp2=xydata
   miss=where(tmp2 lt minval)
;   if miss(0) ne -1 then tmp2(miss)=lev0
   contour,tmp2,xarr,latitude,/overplot,min_value=minval $
       ,levels=clevs,/cell_fill
endif 
if lfill ne 1 then begin
   contour,xydata,xarr,latitude,/overplot,min_value=minval $
       ,levels=clevs,c_linestyle=lines,c_labels=labs,c_thick=cthick
;       ,levels=[0.],c_labels=[1],c_thick=cthick
endif
if projchoice ne 3 then map_continents
map_grid,lons=along,lats=alatg
; map_grid,lons=along,lats=alatg,color=255

return 
end
; *********************************************************
pro monocheck,xa,ya,mu,yin,ivarlab
;
;     Search for undershoots, overshoots or flat areas in the
;     interpolated function, YIN, on the latitude grid, MU.
;     If they occur between neighbouring points, XA, discretised by
;     PV/theta, YA, then apply linear interpolation over this entire
;     interval.
;
;     NOTE: U is not monotonic and naturally has a maximum at the jet.
;     Linear interpolation is not applied over the interval containing
;     the maximum YLAT.
;
;     NOTE: XA,YA ordered from south to north, but MU, YIN
;     start at North Pole, so arrays must be reversed.
;
n=n_elements(xa)
nlat=n_elements(mu)
xlat=reverse(mu)
ylat=reverse(yin)
ymax=max(ylat,jmax)
jlo=0
jhi=1
;
;     Find regular latitude points falling within same interval.
;     Check for undershoot, overshoot or flat section.
;
for k=0,n-2 do begin
    iend=0
    iover=0
    dx=xlat-xa(k)
    if dx(0) lt 0. then begin
        jlo=where(dx ge 0.)
        jlo=jlo(0)
        jhi=jlo
    endif
    if jlo eq -1 then begin
        iend=1
    endif
    while iend eq 0 and jhi le nlat-1 do begin
        dx=xlat(jhi)-xa(k)
        ex=xlat(jhi)-xa(k+1)
        dy=ya(k+1)-ya(k)
        ey=ylat(jhi)-ylat(jhi-1)
        if dx*ex gt 0. then begin
            iend=1
        endif else begin
            if dy*ey le 0. then begin
                iover=1
            endif
        endelse
        jhi=jhi+1
    endwhile
    jhi=jhi-1
;
;     Do not apply adjustment over interval nearest equator
;     if input switch IVARLAB > 0.
;
    if ivarlab gt 0 then begin
        if jlo eq 0 then begin
            iover=0
        end
    endif
;
;     Do not apply adjustment over interval containing maximum YLAT,
;     provided that it does not occur in the interval nearest the
;     ground. Used for zonal flow only.
;
    if ivarlab eq 2 then begin
        if (((JHI-JMAX)*(JLO-JMAX) LT 0) AND (K GT 0)) then begin
            if iover eq 1 then begin
;
;     Chop off overshoots rather than using linear interpolation.
;
                ytest=ylat(jlo:jhi-1)
                ymax=max([ya(k),ya(k+1)],min=ymin)
                jover=where(ytest gt ymax)
                if jover(0) ne -1 then begin
                    ytest(jover)=ymax
                endif
                jover=where(ytest lt ymin)
                if jover(0) ne -1 then begin
                    ytest(jover)=ymin
                endif
                ylat(jlo:jhi-1)=ytest
            endif
            iover=0
        endif
    endif
;
;     If overshoot has occurred then replace interval with
;     linear interpolation.
;
    if iover eq 1 and jhi gt jlo then begin
        if jlo eq 0 then begin
;
;     Fit curve that tends to zero slope at equator.
;
            coeff=(ya(k+1)-ya(k))/(xa(k+1)^4)
            ylat(jlo:jhi-1)=ya(k)+coeff*(xlat(jlo:jhi-1)^4)
        endif else begin
            slope=(ya(k+1)-ya(k))/(xa(k+1)-xa(k))
            ylat(jlo:jhi-1)=ya(k)+slope*(xlat(jlo:jhi-1)-xa(k))
        endelse
    endif
    jlo=jhi
endfor
yin=reverse(ylat)
        
return
end
; *********************************************************
pro intj2ths,thlev,thetas,fj,feq,thsmin,thsmax,emus,f
;
; Interpolate surface field from latitude grid (ordered from NP) 
; to theta levels using theta_s distribution.
;
; print,thetas
nthlev=n_elements(thlev)
nlat=n_elements(thetas)
;
; Start at lowest thlev (nearest north pole)
;
for m=0,nthlev-1 do begin
    thm=thlev(m)
;
; Find indices just north and south of theta=theta_m
;
    js=where(thetas gt thm)
    js=js(0)
    jn=js-1
    if jn ge 0 then begin
        rk=(thm-thetas(jn))/(thetas(js)-thetas(jn))
        f(m)=(1-rk)*fj(jn)+rk*fj(js)
    endif else begin
        if js eq 0 then begin
;
; thm less than lowest thetas value. Assume f=0 at NP.
;
            f(m)=0.
            if (thm gt thsmin) then begin
                rk=(thm-thsmin)/(thetas(0)-thsmin)
                f(m)=rk*fj(0)
            endif
        endif else begin
;
; thm greater than highest thetas value. Assume hemispheric
; and use theta_s=thetas_max and f=feq at equator.
;
;            f(m)=feq
            if (thm lt thsmax) then begin
                rk=(thm-thetas(nlat-1))/(thsmax-thetas(nlat-1))
                f(m)=(1-rk)*fj(nlat-1)+rk*feq
            endif
        endelse
    endelse
endfor
; print,f

return
end
; *********************************************************
pro intj2pv,m,thm,trlev,pvm,fj,qmaxth,qminth,emu1,emusurf,f,fsurf,fnp
;
; Interpolate variable f from latitude grid to PV levels
; using PV distribution on isentropic surfaces.
;
ntrlev=n_elements(trlev)
nlat=n_elements(fj)
f(*,m)=fsurf
kran=where(emu1 lt 1 and emu1 gt emusurf)
nk=n_elements(kran)
kmin=kran(0)
kmax=kran(nk-1)
qedge=qminth(m)
if qedge ge 0. then begin
    if emusurf eq 0 then begin
        qedge=0.
    endif else begin
        qedge=trlev(kmin)
    endelse
endif
;
; Start at highest PV (at NP).
;
for k=kmax,kmin,-1 do begin
    qk=trlev(k)
;
; Find latitude indices just north and south of q=Q_k
; Note j runs from NP.
;
    js=where(pvm lt qk and pvm gt qedge)
    js=js(0)
    jn=js-1
    if jn ge 0 then begin
        rk=(qk-pvm(js))/(pvm(jn)-pvm(js))
        f(k,m)=(1-rk)*fj(js)+rk*fj(jn)
    endif else begin
        if js eq 0 then begin
;
; qk greater than all PV values on this isentropic surface
;
            f(k,m)=fnp
; created problems for large qmaxth near the ground so have omitted.
;            if qk lt qmaxth(m) then begin
;                rk=(qk-pvm(0))/(qmaxth(m)-pvm(0)) 
;                f(k,m)=(1-rk)*fj(0)+rk*fnp
;            endif
        endif else begin
;
; qk less than all PV values on surface
;            
            if qk gt qedge then begin
;                print,'qk > qedge ',m,thm,qedge,k
;                rk=(qk-qedge)/(pvm(nlat-1)-qedge)
;                f(k,m)=(1-rk)*fsurf+rk*fj(nlat-1)
                f(k,m)=fsurf
            endif
        endelse
    endelse
endfor
if kmin gt 0 then begin
    f(0:kmin-1,m)=fsurf
endif
if kmax lt ntrlev-1 then begin
    f(kmax+1:ntrlev-1,m)=fnp
endif

return
end
; *********************************************************
pro filter121,nlat,nthlev,vin,jlower
;
; Apply 1-2-1 filter in the vertical to variable var.
;
var=vin
alpha=0.25
beta=1.-2*alpha
for m=1,nthlev-2 do begin
    jl=jlower(m-1)
    var(0:jl,m)=beta*vin(0:jl,m)+alpha*(vin(0:jl,m-1)+vin(0:jl,m+1))
endfor
vin=var

return
end
; *********************************************************
function trap_uneven,x,f
;
; Use trapezoidal rule to integrate function f known at
; irregular coordinates x.
;
n=n_elements(x)
sum=0.
for i=0,n-2 do begin
   sum=sum+0.5*(f(i)+f(i+1))*(x(i+1)-x(i))
endfor 

return,sum
end


; *********************************************************
; *********************************************************
; Main program backwa_pvinv
;
; *********************************************************
; *********************************************************
hard=0

; Operating system
case (!version.os_family) of
    'Windows': lcomputer=1
    'Unix': lcomputer=0
    else: lcomputer=0
 endcase

; Running by sbatch? Modifications added by Hannah Croad
sbatch_running = (getenv('SLURM_JOB_ID') ne '')
if sbatch_running then lcomputer=2

lcomputer = 2

if lcomputer eq 0 then begin
   set_plot, 'x'
   device,get_visual_name=visname
   device,get_visual_depth=visdepth
   device,get_decomposed=decom
   print,'Display type = ',visname,visdepth,'   decomposed = ',decom
   case visname of
      'PseudoColor': device,pseudo_color=visdepth
      'TrueColor': device,true_color=visdepth,decomposed=0
      else: device,pseudo_color=visdepth
   endcase
   print,'Selected UNIX workstation.'

endif else if lcomputer eq 1 then begin
   set_plot, 'win'
   device, decomposed=0
   print,'Selected Microsoft Windows.'

endif else if lcomputer eq 2 then begin
   set_plot, 'z'
   print, 'Running in batch processing mode (sbatch)'
endif

; ctab=0
ctab=38
loadct,ctab ; Rainbow18
; loadct,0    ; Greyscale
; loadct,4    ; Different rainbow (black bg, yellow text)
tvlct, red, green, blue, /get
r0 = red[0] & g0 = green[0] & b0 = blue[0]
if hard ne 1 and r0 eq 0 and g0 eq 0 and b0 eq 0 then begin
   red[0] = 255B & green[0] = 255B & blue[0] = 255B
   red[255]=0B & green[255]=0B & blue[255]=0B
   tvlct, red, green, blue
endif

;
; Report on the sbatch mode operation
;
if sbatch_running then begin
   if N_ELEMENTS(COMMAND_LINE_ARGS()) lt 11 then begin
      print, 'Usage: idl -e ".run ''backwa_pvinv''" -args <arg1_datasource> <arg2_mode> <arg3_wavenumber> <arg4_startdate> <arg5_nhist> <arg6_length> <arg7_mode1input> <arg8_mode1output> <arg9_mode2input> <arg10_mode2output> <arg11_nthlev>'
      return
   endif
   args = COMMAND_LINE_ARGS()
endif else begin
   ;
   ; See code below for the meaning of each input parameter
   ;
   ldatasrc = 0
   lmode = 1
   mselect = 0
   DATE_START='2010012218'
   FREQ_DAYS=0.25
   nhist=round(FREQ_DAYS*24)
   N_DAYS=0
   tplus=round(N_DAYS*24)
   nthlev_strat=43 ;levels even in pseudoheight - total of 131 theta levels
   
   args = strarr(11)
   args[0] = string(ldatasrc)
   args[1] = string(lmode)
   args[2] = string(mselect)
   args[3] = DATE_START
   args[4] = string(nhist)
   args[5] = string(tplus)
   args[10] = string(nthlev_strat)
endelse

;
; Ask where the data we're reading comes from
;
ldatasrc = LONG(args[0])
; Print out selected data source for verification
case ldatasrc of
   0: print, 'Data source: ERA-I selected'
   1: print, 'Data source: iGCM selected'
   else: print, 'Invalid data source selected'
endcase

;
; Select the mode
;
lmode = LONG(args[1])
; Print out selected mode for verification
case lmode of
   1: print, 'Mode: Calculate integrals'
   2: print, 'Mode: Calculate wave activity given MLM background state'
   3: print, 'Mode: Calculate wave activity relative to neutral state'
   4: print, 'Mode: Calculate integrals, spawn PV inverter and calculate wave activity'
endcase
smode=''
if lmode eq 3 then begin
    smode='N'
endif

lbsevol=1 ; 0 - use Day 0 background, 1 - use current MLM background
lnhonly=0 ; 0 - calculate BS and WA globally, 1 - NH only
ldiscret=2 ; choose level discretisation used for theta and PV
           ; 0 - include Underworld
           ; 1 - Middle and Overworld only
           ; 2 - levels regular in pseudo-height and Lait PV
           ; 3 - appropriate PV-spacing for LC expts
lpvplot=0  ; 0 - no PV map, 1 - PV or PV anomaly, 2 - difference in PV calcs
           ; 3 - also plot pressure at boundaries and PV near model top
           ; 4 - map of pseudomomentum density
           ; 5,6 - map of u,v (full field after any filtering)
           ; 7  - map of lower boundary theta and pressure/streamfunction
           ; 11 - map of PV and longitude of max(PV) versus theta_m (at 65N)
           ; 12 - map of wave-relative streamfunction on theta

ldterms=1 ; 0 - include d-terms in interior WA, 1 - separate d-terms 
mselect=LONG(args[2]) ; 0 - do not filter perturbations to one zonal wavenumber
                      ; >0 - filter to obtain monochromatic perturbation in eta-coords
lffteta=1 ; 1 - FFT fields in eta before interpolation to theta; 0 - do not

c2fix=1.
; cphase=11.888135 ; NM phase speed for LC1 (deg/day)
; cphase=5.630530 ; NM phase speed for LC2 (deg/day)
cphase=0. ; stationary reference frame for streamfunction

;
; Select start date
;
; indate=lonarr(1)
indate = LONG(args[3])
print, 'Start date of run: ', indate
year = STRMID(indate, 2, 4)

; read,indate,prompt='Enter start date of run (YYYYMMDDHH) '
istdate=indate(0)
iday=(istdate-10000*(istdate/10000))/100

expstem='di'
; expstem='lc1'
fstemhi=expstem+'head'
fstemoi=expstem+'out'
if lbsevol eq 0 then begin
    fstembi=fstemoi+'Z'
endif else begin
    fstembi=fstemoi
endelse
expid='_NHANMW'+smode+'_'
; expid='_LC1'+smode+'_'
; expid='_IGCM'+smode+'_'
scrap=''

;
; Specify length
;
tplus = LONG(args[5])
print, 'Length of run is', tplus, ' hours'
tplus=round(tplus)
nhist=LONG(args[4])
print, 'Interval between records is', nhist, ' hours'
npts=tplus/nhist+1
print, 'Number of time points in requested time series is ',npts
twodt=nhist*2./24.

idatadate=istdate   ; was: idatadate=add2date(istdate,tplus)
print,'Calculating background state until ',idatadate

lheating=0
if lheating eq 1 then begin
    print,'WARNING: for heating calculation need records at'
    iplusdate=istdate ; was: iplusdate=add2date(idatadate,nhist)
    print,istdate,' and ',iplusdate
    print,' '
endif

;
; Directories
;
if sbatch_running then begin
   fpathinp=args[6]+'/'
   fpathbs=args[7]+'/'
   fpathpv=args[8]+'/'
   fpatho=args[9]+'/'
endif else begin
   fpathstem='/storage/research/s2senm/df174909/ENM_data/'
   ; File path of input data for mode 1
   ; fpathinp=fpathstem+'step1_convert_analyses/'+year+'/'
   fpathinp='/storage/research/diamet/swrmethn/inv3/tradv/diag_2009/'
   ; File path of output data for mode 1
   ; fpathbs=fpathstem+'step2_integrals/'+year+'/'
   fpathbs=fpathinp+'wadiag/'
   ; File path of input data for mode 2
   fpathpv=fpathstem+'step3_new/'+year+'/'
   ; fpathpv=fpathinp+'pvinv/'
   ; File path of output data for mode 2
   fpatho=fpathstem+'step4_waden/'+year+'/'
endelse
fpathi=fpathinp
fpathjob=fpathi

;
; Set constants.
;
radea=6371299.
omega=7.292e-5
ga=9.80665
rdgas=287.
kappa=2./7.
cp=rdgas/kappa
expo=(1.-kappa)/kappa
rho0=1.2   ; lower boundary density used for plotting (kg m^-3)
p00=1.e5   ; used in definition of potential temperature 
pvu=1.e6   ; used to convert PV from SI units to PVU
thref=380. ; used in definition of modified (Lait) PV (value nr tropical tpp)
           ; also used to define strat for some modifications (like PV chop).
           ; Does not affect lide cycles.
tref=205.  ; used for isothermal atm in conversion to height ( " )
thvalplot=320. ; used for plotting PV anomalies
lbound=1   ; 0 - lower boundary at eta=1
           ; 1 - at eta-level nbound from ground.
           ; 11 - lower boundary given by isentropic level mlb 
nbound=1   ; >=1  (used for lbound=1)
mlb=2      ; >=1  (used for lbound>10)
lmovelb=0  ; 1 - shift lower boundary of 3D state down by orographic height.
           ;     extrapolate variables into new volume down to the geoid.
           ;     Introduced to avoid re-arrangement of hot spots in theta_s.
           ; 0 - no shift.
xtropic=0.1 ; only used here in post-processing of wave activity (mu < xtropic)
case ldiscret of
0: mlb=18
1: mlb=2
2: mlb=1
3: mlb=1
endcase
;
; Parameters associated with removing high PV at low levels, high latitudes 
;
emucut=0.985 ; size of PV spikes to cut off next to pole (all levels)
areacut=0.5*(1-emucut)
case ldatasrc of
   0: begin
      mumod=0.760 ; defines polar LT where isentropic layers have emus > mumod.
                                ; modification occurs in this region
                                ; where d(PV)/d(theta) < 0.
   end
   1: begin
      mumod=1.0                 ; value used for life cycles and IGCM expts
   end
endcase
mband=2      ; number of theta levels to spread lowest level density,
             ; used in part where column mass is corrected to match psurf.
;
; The following parameters only used in polar LT where theta < thetas(mumod)
;                                                          m < mmed
pvtrop=0.75  ; value (~f/r) to chop tropospheric PV to avoid diabatic spike
pvtrop=pvtrop/pvu
areachop=2*areacut
emuchop=1-2*areachop  ; defines polar LT corner where PV is limited < pvtrop
fthres=1.0*2*omega    ; limits on magnitude of relative vorticity
;
; Latitude band for outputting meridional wind for Hovmoellers
;
phibands=45.
phibandn=60.

atmmass=4.*!dpi*radea*radea*p00/ga
afac=1./(4.*!dpi)
gasfac=p00/(rdgas*ga)
nclev=32
massinc=4.
clevs=findgen(nclev)*massinc
clines=clevs
clines(*)=0
clines(0)=1
clevs(0)=1.e-4
tiny=1.e-3
;
; Selected isentropic levels for background state diagnostics.
;
if ldatasrc eq 1 then begin
   ;
   ; Setup for IGCM data
   ;
   moct=1
   thlev0=244.
   thmid=300.
   thtoplev=400.
   dthunder=2.
   dth=2.
   nunder=round((thmid-thlev0)/dthunder)
   thunder=thlev0+findgen(nunder)*dthunder
   nover=round((thtoplev-thmid)/dth)+1
   thover=thmid+findgen(nover)*dth
   thlev=[thunder,thover]
   nthlev=n_elements(thlev)
endif else begin
   case expstem of 
      'th': begin
         moct=1
         thlev0=330.
         thlev=[330.,350.,380.,425.,480.,550.,620.,700.,800.,900.,1000.,1150.$
      ,1300.,1500.,1700.,1900.,2200.,2600.]
   end
      'di': begin
         moct=1
         case ldiscret of
            0: begin
               thlev0=240.
               thmid=376.
               thmid2=thmid
               thtoplev=540.
               dthunder=4.
               dthtop=4.
               thvalplot=320.
     ; Uniform dtheta across Underworld to thmid.
               nunder=round((thmid-thlev0)/dthunder)
               thunder=thlev0+findgen(nunder)*dthunder
               nover=round((thtoplev-thmid2)/dthtop)+1
               thover=[thmid2+findgen(nover)*dthtop]
               thlev=[thunder,thover]
            end
            1: begin
               thlev0=296.
               thmid=376.
               thmid2=700.
               thtoplev=1488.
               dthunder=4.
               dthtop=32.
               thvalplot=312.
     ; Uniform dtheta across Underworld to thmid.
               nunder=round((thmid-thlev0)/dthunder)
               thunder=thlev0+findgen(nunder)*dthunder
     ; Linearly increasing dtheta from thmid to thtop.
               slope=(dthtop-dthunder)/(thmid2-thmid)
               thval=thmid
               thover=thmid
               while thval+dthtop lt thmid2 do begin
                  dth=(thval-thmid)*slope+dthunder
                  thval=thval+dth
                  thover=[thover,thval]
               endwhile
               nover=round((thtoplev-thmid2)/dthtop)+1
               thover=[thover,thmid2+findgen(nover)*dthtop]
               thlev=[thunder,thover]
            end
            2: begin
               thlev0=218.
;       thlev0=240.5
;       thmid=400.
;       dthunder=4.
;       nunder=round((thref-thlev0)/dthunder)
               thmid=450.
               thlow=320.
               dthunder=1.5
;
; Top levels (above thmid) are evenly spaced in pseudo-height.
;
;       dz=1.0 & nthlev=43 ; as used for 2009/2010  (131 levels)
;       dz=1.0 & nthlev=41 ; as used for jj2007     (129 levels)
;       dz=1.0 & nthlev=30 ; lower theta-top (118 levels)
               dz=1.0 & nthlev=LONG(args[10])
               hden=6.5         ; density scale height in km (RT/g for T=222K)
               zscal=hden/kappa ; exponential height scale for theta
               z0=zscal*alog(thmid/thref)
               zarr=findgen(nthlev)*dz+dz+z0 ; to go with 1 or 2K spacing
               thlev=thref*exp(zarr/zscal)
;
; Use regular d(theta) in the Underworld and then expand mid-range
; gaps to blend in which pseudo-height levels in Overworld.
;       
               nunder=round((thlow-thlev0)/dthunder)
               thunder=thlev0+findgen(nunder)*dthunder
               dthmid=thlev(1)-thlev(0)
;       bexpand=0.83*(dthmid-dthunder)/(thmid-thlow) ; goes with dthunder=2K
;       bexpand=0.94*(dthmid-dthunder)/(thmid-thlow) ; goes with dthunder=1K
               bexpand=0.982*(dthmid-dthunder)/(thmid-thlow) ; goes with dthunder=1.5K
               th=thlow
               dtcur=dthunder
               while th lt thmid do begin
                  thunder=[thunder,th]
                  dtcur=dthunder+bexpand*(th-thlow)
                  th=th+dtcur
               endwhile
               ;
               ; Append the Overworld levels to lower theta levels
               ;
               if lbound lt 10 then begin
                  if thmid-thref eq 20. then begin
                     thlev=[thunder,thref+1,thref+9,thlev]
                  endif else if dthunder eq 2. then begin
                     thlev=[thunder,thlev]
                  endif else if dthunder eq 1. then begin
                     thlev=[thunder,thlev]
                  endif else begin
                     thlev=[thunder,thlev]
                  endelse
               endif
               dth=thlev(nthlev-1)-thlev(nthlev-2)
            end
         endcase
         nthlev=n_elements(thlev)
      end
      else: begin
         moct=6
         thlev0=244.
         thmid=300.
         thtoplev=400.
         dthunder=2.
         dth=2.
         nunder=round((thmid-thlev0)/dthunder)
         thunder=thlev0+findgen(nunder)*dthunder
         nover=round((thtoplev-thmid)/dth)+1
         thover=thmid+findgen(nover)*dth
         thlev=[thunder,thover]
         nthlev=n_elements(thlev)
      end
   endcase
endelse

nthlev=n_elements(thlev)
thlevrev=reverse(thlev) 
thlevh=dblarr(nthlev+1)
dtharr=dblarr(nthlev)
dthothh=dblarr(nthlev)
thlevh(0)=thlev(0)-(thlev(1)-thlev(0))/2.
for l=1,nthlev-1 do begin
   thlevh(l)=(thlev(l-1)+thlev(l))/2.
endfor
thtop=thlev(nthlev-1)+(thlev(nthlev-1)-thlev(nthlev-2))/2.
thlevh(nthlev)=thtop
for l=0,nthlev-2 do begin
    dtharr(l)=thlev(l+1)-thlev(l)
    dthothh(l)=dtharr(l)/thlevh(l+1)
endfor
l=nthlev-1
dtharr(l)=thtop-thlev(l)
dthothh(l)=dtharr(l)/thlevh(l+1)
; thz=0.001*rdgas*tref*alog(thlev/tref)/(kappa*ga)
hrho=6.5                        ; density height scale in km
htheta=hrho/kappa
thz=htheta*alog(thlev/thref)

lait2pv=findgen(nthlev)
lait2pv(*)=1.
llait=1

print,nthlev,' isentropic levels (K)'
print,thlev

;plot,findgen(nthlev),thz,xstyle=1,xrange=[0,120]
;stop
;
; Selected modified PV levels.
;
case expstem of 
'th': begin
   ntrlevst=25
   tr0=0.
   citr=2.
   trlevst=tr0+findgen(ntrlevst)*citr
end
'di': begin
    case ldiscret of
    0: begin
      ; as used for summer 2007 for troposphere and LS
        ntrlevst=46
        trlevst=5.+1.*findgen(ntrlevst)
        ntrlevtpp=25
        lnqmin=-1.
        lnqmax=alog(trlevst(0))
        dlnq=(lnqmax-lnqmin)/float(ntrlevtpp-1)
        lnq=lnqmin+dlnq*findgen(ntrlevtpp-1)
        trlevtpp=exp(lnq)
        ntrtrop=10
        trlevtrop=(trlevtpp(0)/float(ntrtrop))*findgen(ntrtrop)
        trlevst=[trlevtrop,trlevtpp,trlevst]
    end
    1: begin
      ; as used for winter 2009 to stratopause
        ntrlevtpp=56
        lnqmin=-1.
        lnqmax=10.
        dlnq=(lnqmax-lnqmin)/float(ntrlevtpp-1)
        lnq=lnqmin+dlnq*findgen(ntrlevtpp-1)
        trlevtpp=exp(lnq)
        ntrtrop=10
        trlevtrop=(trlevtpp(0)/float(ntrtrop))*findgen(ntrtrop)
        trlevst=[trlevtrop,trlevtpp]
    end
    2: begin
      ; expecting to use Lait PV 
        ntrlevst=56
        trlevst=5.+1.*findgen(ntrlevst)
        ntrlevtpp=14
        lnqmin=-1.
        lnqmax=alog(trlevst(0))
        dlnq=(lnqmax-lnqmin)/float(ntrlevtpp-1)
        lnq=lnqmin+dlnq*findgen(ntrlevtpp-1)
        trlevtpp=exp(lnq)
        ntrtrop=5
        trlevtrop=(trlevtpp(0)/float(ntrtrop))*findgen(ntrtrop)
        trlevst=[trlevtrop,trlevtpp,trlevst]
        ioverw=where(thlev ge thref)
;        ioverw=where(thlev ge 0)
        lait2pv(ioverw)=(thlev(ioverw)/thref)^(4.5)
    end
    3: begin
       ntrlevst=32
       trlevst=4.5+0.5*findgen(ntrlevst)
       ntrlevtpp=25
       lnqmin=-1.
       lnqmax=alog(trlevst(0))
       dlnq=(lnqmax-lnqmin)/float(ntrlevtpp-1)
       lnq=lnqmin+dlnq*findgen(ntrlevtpp-1)
       trlevtpp=exp(lnq)
       ntrtrop=10
       trlevtrop=(trlevtpp(0)/float(ntrtrop))*findgen(ntrtrop)
       trlevst=[trlevtrop,trlevtpp,trlevst]
    end
    endcase
end
else: begin
   case ldiscret of
   0: begin
       ntrlevst=16
       strat0=5.
       citrstrat=1.
       troptop=1.
       citrtrop=0.05
       ntrtrop=round(troptop/citrtrop)+1
       trlevstrat=strat0+findgen(ntrlevst)*citrstrat
       trlevtrop=findgen(ntrtrop)*citrtrop
       trlevst= $
          [trlevtrop,1.1,1.25,1.5,1.75,2.,2.25,2.5,3.,3.5,4.,4.5,trlevstrat]
   end
   1: begin
       ntrlevst=41
       lnqmin=-1.
       lnqmax=3.
       dlnq=(lnqmax-lnqmin)/float(ntrlevst-1)
       lnq=lnqmin+dlnq*findgen(ntrlevst)
       trlevst=exp(lnq)
       ntrtrop=10
       trlevtrop=(trlevst(0)/float(ntrtrop))*findgen(ntrtrop)
       trlevst=[trlevtrop,trlevst]
   end
   2: begin
       ntrlevst=32
       trlevst=4.5+0.5*findgen(ntrlevst)
       ntrlevtpp=25
       lnqmin=-1.
       lnqmax=alog(trlevst(0))
       dlnq=(lnqmax-lnqmin)/float(ntrlevtpp-1)
       lnq=lnqmin+dlnq*findgen(ntrlevtpp-1)
       trlevtpp=exp(lnq)
       ntrtrop=10
       trlevtrop=(trlevtpp(0)/float(ntrtrop))*findgen(ntrtrop)
       trlevst=[trlevtrop,trlevtpp,trlevst]
   end
   3: begin
       ntrlevst=32
       trlevst=4.5+0.5*findgen(ntrlevst)
       ntrlevtpp=25
       lnqmin=-1.
       lnqmax=alog(trlevst(0))
       dlnq=(lnqmax-lnqmin)/float(ntrlevtpp-1)
       lnq=lnqmin+dlnq*findgen(ntrlevtpp-1)
       trlevtpp=exp(lnq)
       ntrtrop=10
       trlevtrop=(trlevtpp(0)/float(ntrtrop))*findgen(ntrtrop)
       trlevst=[trlevtrop,trlevtpp,trlevst]
   end
   endcase
end
endcase
ntrlevst=n_elements(trlevst)
trlevsth=dblarr(ntrlevst+1)
trlevsth(0)=trlevst(0)-(trlevst(1)-trlevst(0))/2.
for l=1,ntrlevst-1 do begin
    trlevsth(l)=(trlevst(l-1)+trlevst(l))/2.
endfor
trlevsth(ntrlevst)=trlevst(ntrlevst-1) $
  +(trlevst(ntrlevst-1)-trlevst(ntrlevst-2))/2.

print,' '
laittest=where(lait2pv ne 1)
if laittest(0) ne -1 then begin
    print,' Using modified (Lait) definition of PV throughout.'
    llait=1
endif
print,ntrlevst,' modified PV levels (PVU)'
print,trlevst
trlevst=trlevst/pvu
trlevsth=trlevsth/pvu
emus=dblarr(nthlev)
emus(*)=0.
emu=dblarr(ntrlevst,nthlev)
emu(*,*)=0.
;
; Read potential temperature for first time level 
; (if lheating=1 this will be the timestep before "current").
;
read_trdata,istdate,fpathi,fstemhi,fstemoi,0 $
           ,itratypi,longitude,latitude,aeta,beta,psurf,datai

nlon = n_elements(longitude)
nlat = n_elements(latitude)
nlp = n_elements(aeta)
nlev = nlp-1
;
; Define Gaussian grid and midpoints defining staggered grid for winds.
;
longitude=longitude/moct
longr=longitude*!dpi/180.
dlon=2.*!dpi/float(nlon)  ; Note: do not divide by moct for global weighting
dlongrid=dlon/moct
longh=dblarr(nlon+1)
longh(0)=longr(0)-dlongrid/2.
longh(1:*)=longr(*)+dlongrid/2.
;print,longh*180./!dpi
seglon=360./moct
segtest=0.5*seglon*!dpi/180.
mu=sin(latitude*!dpi/180.)
cosj=cos(latitude*!dpi/180.)
muh=dblarr(nlat+1)
muh(0)=1.
for j=1,nlat-1 do begin
   muh(j)=0.5*(mu(j-1)+mu(j))
endfor
;
nlatend=nlat-1
nlatwa=nlat
if mu(nlat-1) lt 0. then begin
; Global data
   nlatnh=nlat/2
   muh(nlat)=-1.
   jeqp=nlatnh-1
   if lnhonly eq 1 then begin
;
;  Determine whether WA calculation is global or not
;
      nlatwa=nlatnh
;
;  Fix to speed up calculation by omitting SH extratropics
;
      print,''
      print,'WARNING: Omitting SH extratropics in calculations'
      nlatend=nlatnh+20
;      nlatend=nlatnh
   endif
endif else begin
; Hemispheric data
   nlatnh=nlat
   muh(nlat)=0.
   jeqp=nlat-1
endelse
dmu=dblarr(nlat)
for j=0,nlat-1 do begin
   dmu(j)=muh(j)-muh(j+1)
endfor
;print,'dmu = ',dmu
;
boxarea=dblarr(nlon,nlat)
for j=0,nlat-1 do begin
    boxarea(*,j)=afac*dlon*dmu(j)
endfor
latwa=latitude(0:nlatwa-1)
muwa=mu(0:nlatwa-1)
muarr=dblarr(nlon,nlatwa)
for i=0,nlon-1 do begin
    muarr(i,*)=muwa
endfor
cfac=!dpi/(180.*86400.)
crel=radea*radea*cphase*cfac*mu(0:nlatwa-1)
print,' system phase speed at equator (m/s) = ',crel(nlatwa-1)
if nlat eq nlatnh then begin
    projchoice=3
endif else begin
    projchoice=1
endelse

jbands=where(latitude lt phibands)
jbands=jbands(0)
jbandn=where(latitude lt phibandn)
jbandn=jbandn(0)
vband=dblarr(nlon,nthlev)
vband(*,*)=0.
print,' '
print,' Latitude band for meridional wind output = ',jbands,jbandn
findth=where(thlev ge thvalplot)
thchoice=findth(0)
;
; Calculate eta values.
;
alevh=aeta+beta
aetaf=dblarr(nlev)
betaf=dblarr(nlev)
daeta=dblarr(nlev)
dbeta=dblarr(nlev)
rden=dblarr(nlev)
rdenh=dblarr(nlev-1)
for l=0,nlev-1 do begin
   aetaf(l)=0.5*(aeta(l)+aeta(l+1))
   betaf(l)=0.5*(beta(l)+beta(l+1))
   daeta(l)=aeta(l+1)-aeta(l)
   dbeta(l)=beta(l+1)-beta(l)
endfor
alev=aetaf+betaf
nisen = n_elements(itratypi)
itest=where(itratypi eq 2)
itheta=itest(0)
itest=where(itratypi eq 3)
ipv=itest(0)      
itest=where(itratypi eq 5)
iu=itest(0)      
itest=where(itratypi eq 6)
iv=itest(0)      
itest=where(itratypi eq 8)
iheat=itest(0)      
thetam=reform(datai[*,*,*,itheta],nlon,nlat,nlev)
;
; Find surface orography (geopotential).
;
lr4out=1
if lr4out eq 1 then begin
    fdata=fltarr(nlon)
    zs=fltarr(nlon,nlat)
endif else begin
    fdata=dblarr(nlon)
    zs=dblarr(nlon,nlat)
endelse
zs(*,*)=0.

; Orography file (fporog+fhead, under another user's home directory) is not
; available here, and is only ever used downstream inside `lmovelb eq 1`
; blocks (lmovelb=0 by default, matching the Python translation, which does
; not implement an orographic lower-boundary shift). So always take the
; already-existing "no orography data" fallback below instead of reading it.
if 0 then begin
endif else begin
    if lmovelb eq 1 then begin
        print,' '
        print,' CANNOT shift lower boundary - assuming zero orography'
    endif
    lmovelb=0
endelse
zsin=zs
zsav=total(zs,1)/nlon
;print,' zsav ='
;print,zsav
zsav(*)=0.

print,' '
case lbound of
1: begin
    print,'Lower boundary for inversion domain set at eta = ',alev(nlev-nbound)
end
11: begin    
    print,'Lower boundary for inversion domain set at theta = ',thlev(mlb) 
end
else: begin
    print,'Lower boundary for inversion domain set at eta = 1'
end
endcase
print,' '
;
; Read potential temperature and advection for current time level.
; (already read if lheating=0)
;
if lheating eq 1 then begin
    idatadate=add2date(istdate,nhist) 
    read_trdata,idatadate,fpathi,fstemhi,fstemoi,0 $
              ,itratypin,lonin,latin,aetain,betain,psurf,datai
endif else begin
    idatadate=istdate
endelse
theta=reform(datai[*,*,*,itheta],nlon,nlat,nlev)
pv=reform(datai[*,*,*,ipv],nlon,nlat,nlev) ; Ertel PV (PVU)
pv=pv/pvu ; converted to SI units
u=reform(datai[*,*,*,iu],nlon,nlat,nlev)
v=reform(datai[*,*,*,iv],nlon,nlat,nlev)
if iheat gt -1 then begin
   thetaadv=reform(datai[*,*,*,iheat],nlon,nlat,nlev)
endif else begin
   thetaadv=theta
   thetaadv(*,*,*)=0.
endelse


;
; Loop over time.
;

print,'CHECKPOINT: start of loop over time'
print, npts

;
; Read theta diagnostics for next record (used in heating calc).
;
   iplusdate=istdate ; previously iplusdate=add2date(idatadate,nhist)
   read_trdata,iplusdate,fpathi,fstemhi,fstemoi,0 $
           ,itratypin,lonin,latin,aetain,betain,psurfp,datai
   zs=zsin ; reset surface geopotential array to value read from data

   if hard eq 1 then begin
       set_plot,'ps'
       device, /landscape, /color, bits_per_pixel=8
       device, file='thlbLC2plot'+strtrim(idatadate,2)+'.ps'
;           loadct,ctab
   endif
   case lpvplot of
       2: begin
           !p.multi=[0,2,1]
           !p.charsize=1.0
       end
       else: begin
           !p.multi=[0,1,1]
           !p.charsize=1.5
           !p.charthick=3.
       end
   endcase
;
; Calculate zonal averages of primary 3D variables.
;
    print,'CHECKPOINT: calculate zonal averages of 3D variables'
   thzav=total(theta,1)/nlon
   uzav=total(u,1)/nlon
   vzav=total(v,1)/nlon
   pszav=total(psurf,1)/nlon
;
; Overwrite 3D fields with zonal averages.
;
   loverzonal=0
   if loverzonal eq 1 then begin
       for l=0,nlev-1 do begin
           for j=0,nlatend do begin
               theta(*,j,l)=thzav(j,l)
               u(*,j,l)=uzav(j,l)
               v(*,j,l)=vzav(j,l)
           endfor
       endfor
       for j=0,nlatend do begin
           psurf(*,j)=pszav(j)
       endfor
   endif

   ;********************************************************************
   if lmode gt 1 then begin

      ; Load background state data
      bs_file = fpathpv+'bs_data_'+string(idatadate, format='(I010)')+'.nc'
      id = NCDF_OPEN(bs_file)

      ; Get dimensions (latitudes, thetas & PV levels)
      ; Read latitudes 
      lats_id = NCDF_VARID(id, 'latitudes')
      NCDF_VARGET, id, lats_id, rlatinv
      ; Read thetas
      thetas_id = NCDF_VARID(id, 'thetas')
      NCDF_VARGET, id, thetas_id, thetas_bs
      ; Read PV levels
      PVs_id = NCDF_VARID(id, 'PVs')
      NCDF_VARGET, id, PVs_id, PVs_bs
      
      ; Read equivalent latitudes output by outer iteration of PV inverter.
      ; equivalent latitudes (PV-theta grid)
      emu_id = NCDF_VARID(id, 'emu')
      NCDF_VARGET, id, emu_id, emu
      ; surface equivalent latitudes (theta grid)
      emus_id = NCDF_VARID(id, 'surf_emu')
      NCDF_VARGET, id, emus_id, emus
      ; Re-order arrays so that theta runs from bottom to top.
      emu=reverse(emu,2)
      emus=reverse(emus)
      ; Ensure monotonicity of equivalent latitudes wrt PV. 
      for l=0,nthlev-1 do begin
          minemu=min(emu(*,l),minel)
          if minel gt 0 then begin
              print,'Forcing monotonic emu on theta = ',thlev(l),minel
              emu(0:minel-1,l)=emus(l)
          endif
      endfor
      ;
      ; Get zonally symmetric background state fields
      ; Read pressure (lat-theta grid)
      pres_id = NCDF_VARID(id, 'p')
      NCDF_VARGET, id, pres_id, prinv
      ; Read isentropic density (lat-theta grid)
      r_id = NCDF_VARID(id, 'r')
      NCDF_VARGET, id, r_id, sigmainv
      ; Read temperature (lat-theta grid)
      T_id = NCDF_VARID(id, 'T')
      NCDF_VARGET, id, T_id, tempinv
      ; Read zonal wind (lat-theta grid)
      u_id = NCDF_VARID(id, 'u')
      NCDF_VARGET, id, u_id, uthinv
      ; Read Ertel PV (lat-theta grid)
      epv_id = NCDF_VARID(id, 'epv')
      NCDF_VARGET, id, epv_id, pvinv
      ; Read surface pressure (lat grid)
      ps_id = NCDF_VARID(id, 'surf_p')
      NCDF_VARGET, id, ps_id, psinv
      ; Read surface temperature (lat grid)
      ts_id = NCDF_VARID(id, 'surf_T')
      NCDF_VARGET, id, ts_id, tsinv
      ; Read surface zonal wind (lat grid)
      us_id = NCDF_VARID(id, 'surf_u')
      NCDF_VARGET, id, us_id, usinv
      ; Read surface geopotential (lat grid)
      zs_id = NCDF_VARID(id, 'surf_gp')
      NCDF_VARGET, id, zs_id, zsinv

      ; Re-order arrays so that theta runs from bottom to top.
      prinv=reverse(prinv,2)
      tempinv=reverse(tempinv,2)
      sigmainv=reverse(sigmainv,2)
      uthinv=reverse(uthinv,2)
      pvinv=reverse(pvinv,2)

      ; Copy arrays read in from MLM statefile obtained
      ; by optimal transport into the bg state arrays used later.
      ; This MLM state has already been interpolated to the ERA latitudes
      ; on theta levels.
      zs0=zsinv
      ps0=psinv
      ts0=tsinv
      us0=usinv
      pr0=prinv
      temp0=tempinv
      sigma0=sigmainv
      uth0=uthinv
      pv0=pvinv
      ;
      ; Need to extend the data to grid points between the critical
      ; latitude used in optimal transport and the pole.
      ; Extend with uniform values.
      ;
      scrit=max(emus)
      alatcrit=asin(scrit)*180./!dpi
      jbeyond=where(latwa gt alatcrit)
      jlast=max(jbeyond)+1
      zs0(jbeyond)=zs0(jlast)
      ps0(jbeyond)=ps0(jlast)
      ts0(jbeyond)=ts0(jlast)
      us0(jbeyond)=us0(jlast)
      for m=0,nthlev-1 do begin
         pr0(jbeyond,m)=pr0(jlast,m)
         temp0(jbeyond,m)=temp0(jlast,m)
         sigma0(jbeyond,m)=sigma0(jlast,m)
         uth0(jbeyond,m)=uth0(jlast,m)
         pv0(jbeyond,m)=pv0(jlast,m)
         smax=max(emu(*,m))
         iover=where(emu(*,m) ge min([smax,scrit]))
         emu(iover,m)=1.
         emu(iover(0),m)=smax
         unexpected=where(sigma0(*,m) eq 0 and pv0(*,m) gt 0)
         if unexpected(0) ne -1 then begin
            print,'Unexpected zero density on theta ',m,thlev(m)
            sigma0(unexpected,m)=sigma0(unexpected(0)-1,m)
         endif
      endfor
      iover=where(emus ge scrit)
      emus(iover)=1.
      ;
      ; Create an array for u on the equator by extending from interior.
      ; Also this array contains u at the lower boundary for
      ; isentropic surfaces that intersect the ground.
      ;
      ebu0=dblarr(nthlev)
      alats=asin(emus)*180./!dpi
      for m=0,nthlev-1 do begin
         below=where(latwa le alats(m))
         if below(0) eq -1 then begin
            ebu0(m)=uth0(nlatwa-1,m)
         endif else begin
            ebu0(m)=uth0(below(0)-1,m)
         endelse
      endfor
      
      ; Interpolate background state defined in theta-coords to
      ; eta-coords using vertical interpolation linear in ln(p).
      thineta=dblarr(nlatwa,nlev)
      thineta(*,*)=0.
      ueta=thineta
      ; Loop over latitude
      for j=0,nlatwa-1 do begin
          pf=aetaf*p00+betaf*ps0(j)
          pth=pr0(j,*)
          ;  Loop over eta-levels
          for l=0,nlev-1 do begin
              pl=pf(l)
              zl=alog(pl)
              m=where(pth lt pl and pth gt 0)
              mtop=m(0)
              if mtop gt 0 then begin
                  mbot=mtop-1
                  if pth(mbot) gt 0 then begin
                      zbot=alog(pth(mbot))
                      ztop=alog(pth(mtop))
                      fac=(zl-zbot)/(ztop-zbot)
                      ofac=1.-fac
                      thineta(j,l)=ofac*thlev(mbot)+fac*thlev(mtop)
                      ueta(j,l)=ofac*uth0(j,mbot)+fac*uth0(j,mtop)
                  endif else begin
                  ; pl is next to ground, using surface in interpolation
                      zbot=alog(ps0(j))
                      ztop=alog(pth(mtop))
                      fac=(zl-zbot)/(ztop-zbot)
                      ofac=1.-fac
                      thineta(j,l)=ofac*ts0(j)+fac*thlev(mtop)
                      ueta(j,l)=ofac*us0(j)+fac*uth0(j,mtop)
                  endelse
              endif else begin
                  if mtop eq -1 then begin
                  ; pl is above top theta-level. Would use top BC in
                  ; interpolation, except p_top has not been output here.
                      mbot=nthlev-1
                      zbot=alog(pth(mbot))
                      thineta(j,l)=thtop
                      ueta(j,l)=uth0(j,mbot)
                  endif else begin
                  ; pl is next to ground, using surface in interpolation
                      zbot=alog(ps0(j))
                      ztop=alog(pth(mtop))
                      fac=(zl-zbot)/(ztop-zbot)
                      ofac=1.-fac
                      thineta(j,l)=ofac*ts0(j)+fac*thlev(mtop)
                      ueta(j,l)=ofac*us0(j)+fac*uth0(j,mtop)
                  endelse
              endelse
          endfor
       endfor

       iuplot=0
       if iuplot eq 1 then begin
           uplot=ueta
           ulevs=-70.+10.*findgen(21)
           thlevs=230.+5.*findgen(41)
           ulabs=ulevs
           ulabs(*)=1
           contour,uplot,latwa,alev $
             ,title='u!i0!n  (m s!e-1!n)' $
             ,xtitle='Equivalent latitude',xrange=[0,90],xstyle=1 $
             ,ytitle='eta',yrange=[1,0],ystyle=1 $
             ,levels=ulevs,/cell_fill
           contour,uplot,latwa,alev $
             ,levels=ulevs,c_labels=ulabs,/overplot
           contour,thineta,latwa,alev $
             ,levels=thlevs,c_labels=ulabs,/overplot
       endif
       ; Background state comes in as u.
       ; Both converted to (u/a)*cos(phi).
       for j=0,nlatwa-1 do begin
          ueta(j,*)=ueta(j,*)*cosj(j)/radea
       endfor
   ; End if lmode > 1
   endif
   ;********************************************************************

;
; Find average of equatorial lowest level temperature over oceans
; assuming that it is ocean where zs < 30m.
; Note that the arrays contain geopotential rather than height.
;
   geothres=ga*30.
   iocean=where(zs(*,jeqp) lt geothres)
   nocean=n_elements(iocean)
   print,nocean,' ocean points'
   thocean=total(theta(iocean,jeqp,nlev-nbound))/nocean
   thoceanmax=max(theta(iocean,jeqp,nlev-nbound))
   print,'equatorial ocean lower boundary theta ',thocean,thoceanmax
;
; Full analysis fields come in as u/(a*cos(phi)).
; Convert units to rad s^-1.
;
   for j=0,nlatend do begin
       u(*,j,*)=u(*,j,*)*(1.-mu(j)^2)/86400.
   endfor
   v=v/86400.
;
   if mselect gt 0 then begin
       print,' '
       print,' Calculate perturbations from background state in eta-coords'
       print,' FFT perturbations at each (lat, eta) to Fourier coeffs'
       print,' Filter to zonal wavenumber mselect ',mselect*moct
       print,' Inverse FFT to monochromatic perturbation in physical space'
       print,' Add back background state to create new 3-D input field'
;
; Calculate perturbations to the background state in eta-coords.
; FFT to Fourier coefficients, select only m=mselect and transform
; back.
; Note that forward FFT transform returns coefficients for
; m=0,1,2,3,...,mg/2-1,-mg/2-1,-mg/2-2,...,-3,-2,-1
; Therefore if only the +ve wavenumber is selected for the inverse
; transform to complex array "result", the physical space field =
; 2*Re{result}.
;
; Note that in this section the zonal average (m=0) is used rather
; than the MLM zonal state.
;
       loverw=lffteta
       if moct gt 1 then begin
;           loverw=0 ; do not over-write eta-level data in LCs
       endif
       pmonly=dcomplexarr(nlon)
       pmonly(*)=0.
       psfour=dcomplexarr(nlatwa)
       psfour(*)=0.
       zsfour=psfour
       tsfour=psfour
       usfour=psfour
       vsfour=psfour
       uifour=dcomplexarr(nlatwa,nthlev)
       uifour(*,*)=0.
       vifour=uifour
       pifour=uifour
       qifour=uifour
       for j=0,nlatwa-1 do begin
           ppert=psurf(*,j) ; NOTE this is ps/p00
;           ppert=psurf(*,j)-ps0(j)/p00 ; NOTE this is ps/p00
           pfour=FFT(ppert,-1)
           psfour(j)=pfour(mselect)
           if loverw eq 1 then begin
               pav=total(ppert)/nlon
               pmonly(mselect)=pfour(mselect)
               pfilt=FFT(pmonly,1)
               psurf(*,j)=2*real_part(pfilt)+pav
;               psurf(*,j)=2*real_part(pfilt)+ps0(j)/p00
           endif
           zpert=zs(*,j)
;           zpert=zs(*,j)-zs0(j)
           zfour=FFT(zpert,-1)
           zsfour(j)=zfour(mselect)
           if loverw eq 1 then begin
               zav=total(zpert)/nlon
               pmonly(mselect)=zfour(mselect)
               pfilt=FFT(pmonly,1)
               zs(*,j)=2*real_part(pfilt)+zav
;               zs(*,j)=2*real_part(pfilt)+zs0(j)
           endif
           l=nlev-nbound
           thpert=theta(*,j,l)
;           thpert=theta(*,j,l)-thineta(j,l)
           tfour=FFT(thpert,-1)
           upert=u(*,j,l)
;           upert=u(*,j,l)-ueta(j,l)
           ufour=FFT(upert,-1)
           vpert=v(*,j,l)
           vfour=FFT(vpert,-1)
;
;      Determine the Fourier coefficients for lower boundary variables
;
           tsfour(j)=tfour(mselect)
           usfour(j)=ufour(mselect)
           vsfour(j)=vfour(mselect)
       endfor
       if loverw eq 1 then begin
          for l=0,nlev-1 do begin
             for j=0,nlatwa-1 do begin
                thtest=thtop-thineta(j,l)
                if abs(thtest) gt 0.2*dtharr(nthlev-1) then begin
                   thpert=theta(*,j,l)
;                   thpert=theta(*,j,l)-thineta(j,l)
                   tfour=FFT(thpert,-1)
                   thav=total(thpert)/nlon
                   pmonly(mselect)=tfour(mselect)
                   pfilt=FFT(pmonly,1)
                   theta(*,j,l)=2*real_part(pfilt)+thav
;                       theta(*,j,l)=2*real_part(pfilt)+thineta(j,l)
                   
                   upert=u(*,j,l)
;                   upert=u(*,j,l)-ueta(j,l)
                   ufour=FFT(upert,-1)
                   uav=total(upert)/nlon
                   pmonly(mselect)=ufour(mselect)
                   pfilt=FFT(pmonly,1)
                   u(*,j,l)=2*real_part(pfilt)+uav
;                       u(*,j,l)=2*real_part(pfilt)+ueta(j,l)
                   
                   vpert=v(*,j,l)
                   vfour=FFT(vpert,-1)
                   vav=total(vpert)/nlon
                   pmonly(mselect)=vfour(mselect)
                   pfilt=FFT(pmonly,1)
                   v(*,j,l)=2*real_part(pfilt)+vav
;                       v(*,j,l)=2*real_part(pfilt)
                endif
             endfor
             if nlat gt nlatwa and loverw eq 1 then begin
;              Fourier filter zonal winds on equator.
                j=nlatwa
                upert=u(*,j,l)
                ufour=FFT(upert,-1)
                uav=total(upert)/nlon
                pmonly(mselect)=ufour(mselect)
                pfilt=FFT(pmonly,1)
                u(*,j,l)=2*real_part(pfilt)+uav
             endif               
          endfor
       endif
    endif else begin
       print,'*********************************'
       print,'mselect=0 - no Fourier filtering'
    endelse
;
; Interpolate (u/a)*cos(phi) and theta onto latitude midpoints. 
;
   print,' Horizontal interpolation of winds'
   uh=dblarr(nlon,nlat,nlev)
   uh(*,*,*)=0.             ; boundary condition at North Pole.
   thh=theta  
   for j=1,nlatend do begin
       uh(*,j,*)=(u(*,j-1,*)+u(*,j,*))/2.
       thh(*,j,*)=(theta(*,j-1,*)+theta(*,j,*))/2.
   endfor
;
; Interpolate (v/a)*cos(phi) and theta onto longitude midpoints.
;
   vlon=dblarr(nlon,nlat,nlev)
   thlon=dblarr(nlon,nlat,nlev)
   for i=0,nlon-2 do begin
       vlon(i,*,*)=(v(i,*,*)+v(i+1,*,*))/2.
       thlon(i,*,*)=(theta(i,*,*)+theta(i+1,*,*))/2.
   endfor
   vlon(nlon-1,*,*)=(v(nlon-1,*,*)+v(0,*,*))/2.
   thlon(nlon-1,*,*)=(theta(nlon-1,*,*)+theta(0,*,*))/2.
;
; Interpolate fields from eta to theta coordinates. 
;
   print,' '
   print,' Vertical interpolation from eta to theta'
;
; First find the range of theta-levels that are covered by the input data.
;
   topminth=min(theta(*,*,0))
   print,' Minimum theta on top eta-level of data = ',topminth
   nthlim=where(thlevh gt topminth)
   if nthlim(0) ne -1 then begin
       nthlim=nthlim(0)-1
   endif else begin
       nthlim=nthlev
   endelse
   print,' Top theta boundary for wave activity analysis = ' $
        ,nthlim,thlevh(nthlim)
   print,' Top full-level in theta coordinates = ',thlev(nthlim-1)
   print,' '
;
; Define arrays for isentropic coordinates
;
   tth=dblarr(nlon,nlat,nthlev)
   tth(*,*,*)=0.
   pvthh=tth
   qth=dblarr(nlon,nlat,nthlev)
   qth(*,*,*)=0.
   qmod=qth
   zetath=qth
   zetamod=qth
   denth=qth
   pvth=qth
   dmdlam=qth
   wa3d=qth

   utermarr=qth ; added 5/10/26 to diagnose conversion to Python
   vtermarr=qth ; added 5/10/26 to diagnose conversion to Python

;   wa3dqmt=qth
;   wa3dqm0=qth
;   wa3dct=qth
;   wa3dc0=qth
   mont=qth
   marr=dblarr(nthlev)
   marr(*)=0.
   uth=dblarr(nlon,nlat+1,nthlev)
   uth(*,*,*)=0.
   vth=dblarr(nlon+1,nlat,nthlev)
   vth(*,*,*)=0.
   ptop=dblarr(nlon,nlat)
   mbotarr=lonarr(nlon,nlat)
   ptop(*,*)=0.
   ztop=ptop
   zlb=ptop
   plb=psurf
   thetas=ptop
   thetalb=ptop
   denths=ptop
   zlbav=dblarr(nlatwa)
   plbav=dblarr(nlatwa)
   print,'*** CHECKPOINT 1: about to enter first j-loop, nlatend=',nlatend,' nlon=',nlon
   for j=0,nlatend do begin
       if j eq 0 then print,'*** CHECKPOINT 2: inside first j-loop, j=0 ***'
;       print,' j=',j
       for i=0,nlon-1 do begin
           tharr=reform(theta[i,j,*],nlev)
           thharr=reform(thh[i,j,*],nlev)
           thlonarr=reform(thlon[i,j,*],nlev)
           ps=psurf(i,j) ; NOTE this is ps/p00
           plev=aetaf+betaf*ps
           if lbound gt 0 then begin
;
; Use the lowest model eta-level as the lower boundary.
;
              thsurf=tharr(nlev-1)
              thlb=tharr(nlev-nbound)
              if thsurf gt thlb then begin
;
; Move statically unstable regions near surface to below lower boundary.
;
                  thlb=thsurf
              endif
              thsurf=thlb
           endif else begin
;
; Extrapolate theta linearly in eta to eta=1 to determine surface theta.
;
              thsurf=tharr(nlev-2)+(1.-alev(nlev-2)) $
                     *(tharr(nlev-1)-tharr(nlev-2))/(alev(nlev-1)-alev(nlev-2))
              thlb=thsurf
           endelse
           thetas(i,j)=thsurf
           if lbound gt 10 then begin
;
; Lower boundary of inversion domain will be isentropic level above ground.
;
              dth=thlevh(mlb+1)-thlevh(mlb)
              thlb=thlev(mlb)-0.25*(mu(j)+1)*dth
              thetalb(i,j)=thlb
           endif else begin
              thetalb(i,j)=thlb
           endelse
;
; Find geopotential on eta-levels and then interpolate to
; top full level in theta-coords (assuming everywhere above ground).
;
           tarr=tharr*(plev^kappa)
           zarr=dblarr(nlev)
           zarr(nlev-1)=zs(i,j)+rdgas*tarr(nlev-1)*alog(ps/plev(nlev-1))
           for l=nlev-2,0,-1 do begin
               zarr(l)=zarr(l+1)+rdgas*0.5*(tarr(l+1)+tarr(l)) $
                                          *alog(plev(l+1)/plev(l))
           endfor
           case lbound of
           1: begin
               zlb(i,j)=zarr(nlev-nbound)
           end
           11: begin
;
; Find geopotential on lower isentropic surface. 
;              
              thgr=where(tharr le thetas(i,j))
              lp=thgr(0)
              if lp ge 0 then begin
                 lm=lp-1
                 if lm lt 0 then begin
                    lp=1
                    lm=0
                 endif
                 r=(thetas(i,j)-tharr(lm))/(tharr(lp)-tharr(lm))
                 omr=1.-r
                 zlb(i,j)=omr*zarr(lm)+r*zarr(lp)
              endif
           end
           else: begin
               print,' Need to specify lower boundary geopotential'
               stop
           end
       endcase
;
; Interpolation to full theta levels for latitude midpoints.
; Loop from bottom to top.
;
           for m=0,nthlev-1 do begin
               thm=thlev(m)
               thgr=where(thharr lt thm)
               lp=thgr(0)
               if lp ge 0 then begin
;
; Model (eta) levels ordered from top to ground.
; Model level (lp) just below chosen theta level exists.
;
                   lm=max([lp-1,0])
                   lp=lm+1
;
; Since 11/08/10 assuming that lm and lp same on horizontal half-grids
; for speed up. Appears to have almost no effect.
;
                   r=(thm-thharr(lm))/(thharr(lp)-thharr(lm))
                   omr=1.-r
                   uth(i,j,m)=omr*uh(i,j,lm)+r*uh(i,j,lp)
               endif
; Interpolation to full theta levels for longitude midpoints.
               thgr=where(thlonarr lt thm)
               lp=thgr(0)
               if lp ge 0 then begin
                   lm=max([lp-1,0])
                   lp=lm+1
                   r=(thm-thlonarr(lm))/(thlonarr(lp)-thlonarr(lm))
                   omr=1.-r
                   vth(i+1,j,m)=omr*vlon(i,j,lm)+r*vlon(i,j,lp)
               endif
; Interpolation to full theta levels for Gaussian grid.
               thgr=where(tharr lt thm)
               lp=thgr(0)
               if lp ge 0 then begin
                   lm=max([lp-1,0])
                   lp=lm+1
                   r=(thm-tharr(lm))/(tharr(lp)-tharr(lm))
                   omr=1.-r
                   pm=omr*plev(lm)+r*plev(lp)
                   pm=max([pm,plev(0)])
                   tth(i,j,m)=thm*(pm^kappa)
; Interpolate Ertel PV onto theta level and convert to modified PV.
                   pvth(i,j,m)=(omr*pv(i,j,lm)+r*pv(i,j,lp))/lait2pv(m)
               endif               
           endfor
           if lp gt 0 then begin
               ztop(i,j)=omr*zarr(lm)+r*zarr(lp)
           endif else begin
               print,'WARNING: top theta-level is beyond range of data'
               print,'j,i,lp = ',j,i,lp
               ztop(i,j)=zarr(0)
           endelse
;
; Adjustment of isentropic data associated with shift in the ground 
; to surface geopotential = zonal average from analysis (zsav). 
; Note: z arrays are geopotential (ga*z). Pressure arrays *p00.
; Used only before finding basic state integrals.
;
           plb(i,j)=plev(nlev-nbound)
           if lmovelb eq 1 then begin
               brunt=1.25e-2 ; N to use in extrapolation underground
               th00=300.
               dthdz=(th00/ga)*brunt*brunt/ga
;              Find extrapolated theta location assuming constant N
               thsref=thetalb(i,j)-dthdz*(zs(i,j)-zsav(j)) 
               thgr=where(thlev gt thsref and thlev lt thetalb(i,j))
               if thgr(0) ne -1 then begin
                   zextrap=(thlev-thsref)/dthdz
                   nextrap=n_elements(thgr)
                   zupper=zlb(i,j)
                   pupper=plb(i,j)
                   m=thgr(nextrap-1)
                   mupper=m+1
                   for n=0,nextrap-1 do begin
                       dz=zextrap(m)-zupper
                       thm=thlev(m)
                       dp=-(dz/(rdgas*thm))*pupper*pupper^(-kappa)
                       pupper=pupper+dp
                       tth(i,j,m)=thm*(pupper^kappa)
                       uth(i,j,m)=uth(i,j,mupper)
                       vth(i,j,m)=vth(i,j,mupper)
                       pvth(i,j,m)=pvth(i,j,mupper)
                       zupper=zextrap(m)
                       m=m-1
                   endfor
                   zlb(i,j)=zlb(i,j)-(zs(i,j)-zsav(j))                   
                   dz=zlb(i,j)-zupper
                   thm=thsref
                   dp=-(dz/(rdgas*thm))*pupper*pupper^(-kappa)
                   plb(i,j)=pupper+dp
                   thetalb(i,j)=thsref
               endif
           endif
       endfor             ; end of longitude loop
;
;      Adjust geopotential and pressure on bottom and top boundaries with the
;      zonal mean?
;       
       if j lt nlatwa then begin
          zlbav(j)=total(zlb(*,j))/float(nlon)
          ztbav=total(ztop(*,j))/float(nlon)
          plbav(j)=total(plb(*,j))/float(nlon)
;          plbav(j)=total(psurf(*,j))/float(nlon)
;          zlb(*,j)=zlb(*,j)-zlbav(j)+zs0(j)
;          plb(*,j)=plb(*,j)-plbav(j)+ps0(j)/p00
;          ztop(*,j)=ztbav
       endif
       for i=0,nlon-1 do begin
           marr(*)=0.
           tharr=reform(theta[i,j,*],nlev)
           tarr=tth(i,j,*)
           thlb=thetalb(i,j)
           ps=psurf(i,j) ; NOTE this is ps/p00
           plev=aetaf+betaf*ps
;           
; Find pressure on top boundary.
; Remember tharr ordered from top to ground.
;
           thgr=where(tharr lt thtop)
           lp=thgr(0)
           if lp ge 1 then begin
              lm=lp-1
              r=(thtop-tharr(lm))/(tharr(lp)-tharr(lm))
              omr=1.-r
              ptop(i,j)=omr*plev(lm)+r*plev(lp)
              ttop=thtop*(ptop(i,j)^kappa)
           endif else begin
              if i+j eq 0 then begin
                  print,'    *** thtop above range of eta-level data'
                  print,i,j,thtop,tharr(0)
              endif
              ptop(i,j)=plev(0)
              ttop=thtop*(ptop(i,j)^kappa)
           endelse
;
; Integrate temperature from top full level in theta-coords to ground 
; to find Montgomery potential on full levels.
;           
           tthtop=tarr(nthlev-1)
           zthtop=ztop(i,j)
           mthtop=cp*tthtop+zthtop
;           if i eq 0 then begin
;               print,'zthtop,tthtop =',zthtop,tthtop
;           endif           
           mtest=where(thlev ge thlb)
           mbot=mtest(0)
           mbp=mbot+1
           m=nthlev-1
           marr(m)=mthtop
           for m=nthlev-2,mbot,-1 do begin
               marr(m)=marr(m+1)-cp*0.5*(tarr(m+1)+tarr(m))*dthothh(m)
           endfor
           mont(i,j,*)=marr(*)
;
; Differentiate M to find isentropic density on theta levels.
;
           thm=thlev(mbot)
           thp=thlev(mbp)
           dth=thp-thm
           mm=marr(mbot)
           mp=marr(mbp)
;
; Use quadratic extrapolation of M to surface for first theta-level
; above ground. Assign this value to surface pseudodensity too.
;
           under=dth*(thp*thm-thlb*thlb)
           dmdth=mp*(thm*thm-thlb*thlb) $
                +mm*(thlb*thlb+thp*thp-2*thp*thm)-zlb(i,j)*dth*dth
           dmdth=dmdth/under
           d2mdth2=(mp*thm-mm*thp+zlb(i,j)*dth)*2./under
           denth(i,j,mbot)=-gasfac*d2mdth2*(dmdth/cp)^expo
;
; Now calculate isentropic density for interior theta-levels.
;
           tho=thm
           mo=mm
           for m=mbp,nthlev-2 do begin
               thm=tho
               tho=thp
               thp=thlev(m+1)
               mm=mo
               mo=mp
               mp=marr(m+1)
               dmdth=(mp-mm)/(thp-thm)
               d2mdth2=((mp-mo)/(thp-tho)-(mo-mm)/(tho-thm))*2/(thp-thm)
               denth(i,j,m)=-gasfac*d2mdth2*(dmdth/cp)^expo
           endfor
           m=nthlev-1
           dmdth=cp*tthtop/thlev(m)  ; Using temperature on top theta level.
           d2mdth2=(cp*ttop/thtop-(mp-mo)/(thp-tho))/(thtop-thlevh(m))
           denth(i,j,m)=-gasfac*d2mdth2*(dmdth/cp)^expo
;
; Overwrite density on band of lowest theta-levels above ground 
; (mbot:minband) such that mass of column from lower boundary to top 
; is respected exactly.
;
           if lbound eq 1 then begin
               minband=mbot+mband  ; define band of modification 
               denint=0.
               for m=minband+1,nthlev-2 do begin
                   denint=denint+denth(i,j,m)*(thlevh(m+1)-thlevh(m))
               endfor
               m=nthlev-1
               denint=denint+denth(i,j,m)*(thtop-thlevh(m))
               massbot=(plb(i,j)-ptop(i,j))*p00/ga-denint
               denth(i,j,mbot:minband)=massbot/(thlevh(minband+1)-thlb)
           endif
           denths(i,j)=denth(i,j,mbot)
       endfor
;
; Calculate zonal derivative of Montgomery potential.
;
       dmdlam(0,j,*)=(mont(1,j,*)-mont(nlon-1,j,*))/(2.*dlon)
       for i=1,nlon-2 do begin
           dmdlam(i,j,*)=(mont(i+1,j,*)-mont(i-1,j,*))/(2.*dlon)
       endfor
       dmdlam(nlon-1,j,*)=(mont(0,j,*)-mont(nlon-2,j,*))/(2.*dlon)
   endfor                      ; end of latitude loop
   print,'*** CHECKPOINT 3: first j-loop finished ***'

   boxw=10
   zlbav=smooth(zlbav,boxw)
   plbav=smooth(plbav,boxw)
   latrev=reverse(latwa) ; order Gaussian grid from equator to NP
   zlbav=reverse(zlbav)  ; re-order zonal average of LB geopotential
   if lmode gt 1 then begin
      zs2a=spl_init(latrev,zlbav,/double)
      ztmp=spl_interp(latrev,zlbav,zs2a,rlatinv,/double) ; put on inverter grid
      zsinv=0.5*(ztmp/ga+zsinv) ; over-write LB height by mixing with analysis.
;      zsinv=ztmp/ga ; Only used to determine the boundary intersection psi0b.
   endif
   print,'zlb (0,0,m)  = ',zlb(0,0)/ga
   
   limplot=[0, -180, 90, 180]
   if lpvplot eq 7 then begin
       tlevs=240.+findgen(30)*4.
       tanom=-10.+findgen(30)*1.
       tlabs=tlevs
       tlabs(*)=1
       thplot=thetalb(*,0:nlatwa-1)
;       tlevs=tanom
;       for j=0,nlatwa-1 do begin
;           thplot(*,j)=thplot(*,j)-ts0(j)
;       endfor
       plbeq=plb(0,nlatwa-1)
       stream=dblarr(nlon,nlatwa)
       for j=0,nlatwa-1 do begin
           stream(*,j)=crel(j)+(plb(*,j)-plbeq)*p00/(2*omega*mu(j)*rho0)
       endfor
       maxmont=max(stream,min=minmont)
       print,'streamfunction max, min',maxmont,minmont
       dmont=(maxmont-minmont)/(nclev-2)
       attrchoice=iv
;       units='streamfunction  (m!e2!ns!e-1!n)'
       units='theta  (K)'
       lev0=260 & dlev=2 & llog=0 & lfill=1 & nozer=0
       lsurf=1 & lisen=0 & etachoice=nlev-nbound
       plotxy,idatadate,thplot,longitude,latwa,alev $
         ,lsurf,lisen,thlev,attrchoice,units $
         ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
       lev0=minmont-2*dmont & dlev=dmont & llog=0 & lfill=0 & nozer=0
       lsurf=1 & lisen=0
;       plotxy,idatadate,stream,longitude,latwa,alev $
;         ,lsurf,lisen,thlev,attrchoice,units $
;         ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice

       goto,skip

       print,' Plotting theta on lower boundary'
       map_set, 90, 0, limit=limplot, /stereographic, /isotropic $
         , title='Theta on lower boundary'
;       contour,thetas-thetalb,longitude,latitude $
;              ,/overplot,/cell_fill,levels=tanom
       contour,thplot,longitude,latwa,/overplot,/cell_fill,levels=tlevs
       contour,thplot,longitude,latwa,/overplot,levels=tlevs,c_labels=tlabs
       map_grid,lons=along,lats=alatg
       map_continents
       print,' type .cont to continue'
       stop

;       plevs=800.+findgen(30)*10.
       plevs=900.+findgen(30)*5.
       panom=-500.+findgen(30)*25.
       plabs=plevs
       plabs(*)=1
       print,' Plotting pressure on lower boundary'
       map_set, 90, 0, limit=limplot, /stereographic, /isotropic $
         , title='Pressure on lower boundary'
;       contour,(psurf-plb)*1000,longitude,latitude $
;              ,/overplot,/cell_fill,levels=panom
       contour,plb*1000,longitude,latitude,/overplot,/cell_fill,levels=plevs
       contour,plb*1000,longitude,latitude,/overplot,levels=plevs,c_labels=plabs
       map_grid,lons=along,lats=alatg
       map_continents
       print,' type .cont to continue'
       stop
   endif

   jpick=where(mu lt emuchop)  ; mu ordered from NP southwards
   jpick=jpick(0)
   chopths=total(thetalb(*,jpick))/nlon
   print,'PV values limited by pvtrop near pole below theta = ',chopths
   jpick=where(mu lt mumod)  ; mu ordered from NP southwards
   jpick=jpick(0)
   if jpick ne -1 then begin
       medianths=total(thetalb(*,jpick))/nlon
       mmed=where(thlev gt medianths)
       mmed=mmed(0)
       print,'PV modified where dP/d(theta) < 0 below theta = ',medianths
   endif else begin
       mmed=0
       print,'No density and PV modification in Arctic LT' 
   endelse
   thsmax=max(thetalb(*,0:nlatwa-1),min=thsmin)
   mmin=where(thlev gt thsmin)
   mmin=mmin(0)
   print,' PV modification within domain bounded by:'
   print,jpick,mu(jpick),mmin,thlev(mmin),mmed,thlev(mmed)
;
; Calculate the isentropic vorticity from the winds.
; Assume that U=V=0 below ground.
; Calculate Ertel PV above ground from isentropic vorticity and density.
;
   print,' Isentropic vorticity calculation'
   vth(0,*,*)=vth(nlon,*,*) ; Periodic boundary conditions.
;
; If domain is global the B.C. U=0 at South pole is used.
; If domain is hemispheric, choose between zero vorticity
; or specified equatorial wind at the equator.
;
; Also interpolate u, v from 1/2-grids to full-grid for
; later calculation of energy.
;
   vthf=dblarr(nlon,nlat,nthlev)
   vthf(*,*,*)=0.
   uthf=vthf
   for i=0,nlon-1 do begin
       vthf(i,*,*)=0.5*(vth(i,*,*)+vth(i+1,*,*))
   endfor
   uterm=dblarr(nthlev)
   vterm=uterm
   print,'*** CHECKPOINT 4: about to enter second j-loop (zetath/qth) ***'
   for j=0,nlatend-1 do begin
       if j eq 0 then print,'*** CHECKPOINT 5: inside second j-loop, j=0 ***'
       secsq=1./(1.-mu(j)^2)
       fcor=2.*omega*mu(j)
       jtest=j*(j-nlatend)
       for i=0,nlon-1 do begin
           thsurf=thetalb(i,j)
           mtest=where(thlev ge thsurf)
           mbot=mtest(0)
           mbotarr(i,j)=mbot
           uterm(*)=0.
           vterm(*)=0.
           utest=mtest
           if jtest ne 0 then begin
               utest=where(uth(i,j+1,*)*uth(i,j,*) ne 0.)
           endif
           vtest=where(vth(i+1,j,*)*vth(i,j,*) ne 0.)           
           uthf(i,j,mtest)=0.5*(uth(i,j,mtest)+uth(i,j+1,mtest))
           uterm(utest)=-(uth(i,j+1,utest)-uth(i,j,utest))/(muh(j+1)-muh(j))
           vterm(vtest)=secsq*(vth(i+1,j,vtest)-vth(i,j,vtest))$
                             /(longh(i+1)-longh(i))
           zetath(i,j,mtest)=2.*omega*mu(j)+uterm(mtest)+vterm(mtest)
           
           utermarr(i,j,mtest)=uterm(mtest)
           vtermarr(i,j,mtest)=vterm(mtest)
           
           qth(i,j,mtest)=zetath(i,j,mtest)/(lait2pv(mtest)*denth(i,j,mtest))
           if mbot gt 0 then begin
               vthf(i,j,0:mbot-1)=0.
               dmdlam(i,j,0:mbot-1)=0.
           endif
;
; In extratropics where lower boundary theta is low
; modify density below the first minimum q searching up from ground.
; Boost isentropic density in the bottom boxes that remain, such that
; column integral of mass is consistent with lower and upper pressure.
; The modified boxes are denoted minel:minband and contain equal density.
; Recalculate PV in this box assuming that vorticity is unchanged.
;
; NOTE: not applied to life cycles since nlat=nlatnh
;
           minel=mbot
           if mbot+1 lt mmed and nlat gt nlatnh then begin
;
;  Modify relative vorticity so that it cannot exceed fthres in
;  Arctic Underworld.
;
               for m=mbot,mmed do begin
                   relvort=zetath(i,j,m)-fcor
                   if abs(relvort) gt fthres then begin
                       signvort=relvort/abs(relvort)
                       vortmod=fcor+fthres*signvort
                       zetath(i,j,m)=vortmod
                       qth(i,j,m)=vortmod/(lait2pv(m)*denth(i,j,m))
                   endif
               endfor
               qbot=qth(i,j,mbot)
               m=mbot+1
               qm=qth(i,j,m)
               while qm le qbot do begin
                   m=m+1
                   qbot=qm
                   qm=qth(i,j,m)
               endwhile
               minel=m-1
               if minel gt mbot then begin
                   minband=minel+mband ; define band of modification 
                   thlb=thsurf
;                   thlb=thlev(minel)-0.5*dtharr(minel)
;                   thetalb(i,j)=thlb
                   denint=0.
                   for m=minband+1,nthlev-2 do begin
                      denint=denint+denth(i,j,m)*(thlevh(m+1)-thlevh(m))
                   endfor
                   m=nthlev-1
                   denint=denint+denth(i,j,m)*(thtop-thlevh(m))
                   massbot=(plb(i,j)-ptop(i,j))*p00/ga-denint
                   minel=mbot
                   denth(i,j,minel:minband)=massbot/(thlevh(minband+1)-thlb)
               endif
           endif
           zetamod(i,j,minel:*)=zetath(i,j,minel:*)
           qmod(i,j,minel:*)=zetath(i,j,minel:*) $
                            /(lait2pv(minel:*)*denth(i,j,minel:*))
       endfor
   endfor
   qth=qmod
   print,'*** ABOUT TO SAVE AT LINE 2392 ***'
   save, filename='/home/users/cq934523/ENMs/empirical-normal-modes/tests/interpolation/idl_output_2010022218.sav', $
      uth, vth, pm, tth, utermarr,vtermarr,$
      qth, zetath, zetamod, denth, $
      thetalb, plb, ptop, mbotarr, thlev, thlevh, mmed, medianths, $
      longitude, latitude, nlatend, nlatnh
   print,'*** SAVE COMPLETED ***'
    stop
;   
;  Define arrays for integral diagnostics.
;
   if ipt eq 1 then begin
      bsmass=dblarr(ntrlevst,nthlev)
      bsmodpv=dblarr(nlatwa,nthlev)
      psmass=dblarr(nthlev)
   endif
   bsmass(*,*)=0.
   bscirc=bsmass
   bsplan=bsmass
   bsarea=bsmass
   bsden=bsmass
   dmudq=bsmass
   bsucosk=bsmass
   uthk=bsmass
   sigmak=bsmass
   pcask=bsmass
   qlength=bsmass
   cdivl=bsmass
   bsmodpv(*,*)=0.
   waden=bsmodpv
   wadenlin=waden
   wagrav=waden
   wad=waden
   wae=waden
   waetmp=wae
   peetmp=wae
   peden=waden
   pew=waden
   pegrav=waden
   peape=waden
   ped=waden
   pee=waden
   fluxphi=waden
   fluxth=waden
   fluxdiv=waden
   qbaryj=bsmodpv
   grlatj=bsmodpv
   psmass(*)=0.
   psarea=psmass
   pscirc=psmass
   psvmom2=psmass
   psvmom3=psmass
   psvmom4=psmass
   qmaxth=psmass
   qminth=psmass
   usm=psmass
   wab=psmass
   peb=psmass
   ueq=dblarr(nthlev)
   ueq(*)=0.
   pzon=dblarr(nlatwa)
   pzon(*)=0.
   zs0=pzon ; over-writing previous value of zs0 read in.
   pet=pzon
   petop=pet
;
;  Find Eulerian zonal average zonal u*cos(phi).
;   
   uthf=radea*uthf
   vthf=radea*vthf
   zavucosj=total(uthf(*,0:nlatwa-1,*),1)/float(nlon)
;
;  Calculate zonal average of pressure on top boundary. 
;  Not necessarily a good boundary condition for background state
;  because a strong polar vortex displaced off the pole can 
;  introduce spurious waves are inconsistent with the monotonic 
;  PV distribution of the re-arranged state.
;
;  So identify polar vortex centre with max(column PV) and rotate
;  coordinates so that pole passes through vortex centre.
;  Calculate zonal average of transformed field without re-gridding
;  data since this is very time consuming.
;
   ist=floor(nthlev*3/4)
   ien=nthlev-1
   qvert=total(qth(*,*,ist:ien),3)
   utop=uthf(*,*,ien)
   test=where(~finite(qvert))
   if test(0) ne -1 then begin
       qvert(test)=0.
   endif
   maxq=max(qvert,kmaxq)
   jmaxq=floor(kmaxq/nlon)
   imaxq=kmaxq-jmaxq*nlon
   lonc=longitude(imaxq)
   latc=latitude(jmaxq)
   print,' Max PV on top theta-level = ' $
        ,maxq,qvert(imaxq,jmaxq),lonc,latc
   lon2d=dblarr(nlon,nlat)
   for j=0,nlat-1 do begin
      lon2d(*,j)=longitude*!dpi/180.
   endfor
   lat2d=dblarr(nlon,nlat)
   for i=0,nlon-1 do begin
      lat2d(i,*)=latitude*!dpi/180.
   endfor
   limplot=[0, -180, 90, 180]
   DDEG=15.
   NLONG=360./DDEG
   NLATG=180./DDEG
   ALONG=FINDGEN(NLONG+1)*DDEG
   ALATG=FINDGEN(NLATG)*DDEG-89.999
   if nlat eq nlatnh then begin
       projchoice=3
   endif else begin
       projchoice=1
   endelse
   if lpvplot eq 3 then begin
       plevs=findgen(20)*0.00001
       plabs=plevs
       plabs(*)=1
       map_set, 90, 0, limit=limplot, /stereographic, /isotropic $
         , title='Interpolated pressure on top boundary'
;       contour,qvert,longitude,latitude,/overplot,/cell_fill
       contour,ptop,longitude,latitude,/overplot,levels=plevs,c_labels=plabs
       map_grid,lons=along,lats=alatg
       stop
   endif
   if lpvplot eq 11 then begin
       qlevs=(findgen(30)-2)*1.
;       qlevs=(findgen(30)+1)*20.
       qlabs=qlevs
       qlabs(*)=1
;       m=113 ; 1411K
       m=118 ; 1758K
       map_set, 90, 0, limit=limplot, /stereographic, /isotropic $
         , title='Ertel PV modified on theta = '+strtrim(round(thlev(m)),2)+'K'
       contour,pvu*qmod(*,*,m),longitude,latitude $
              ,/overplot,/cell_fill,levels=qlevs
;       contour,denth(*,*,m),longitude,latitude $
;              ,/overplot,/cell_fill,levels=qlevs
;       contour,denth(*,*,m),longitude,latitude,/overplot $
;              ,levels=qlevs,c_labels=qlabs
       map_grid,lons=along,lats=alatg
       map_continents
       stop
;
       lonpvmax=fltarr(nthlev)
       jtest=where(latitude gt 60 and latitude lt 75)
       for m=0,nthlev-1 do begin
           qband=total(qmod(*,jtest,m),2)
           maxpv=max(qband,imax)
           lonpvmax(m)=longitude(imax)
       endfor
       iwrap=where(lonpvmax gt 180.)
       lonpvmax(iwrap)=lonpvmax(iwrap)-360.
       plot,lonpvmax,thz $
           ,xtitle='Longitude',ytitle='z (km)'
       stop
   endif
   if lpvplot eq 12 then begin
       monteq=mont(0,nlatwa-1,thchoice)
       if monteq eq 0 then begin
           nonz=where(mont(0,*,thchoice) gt 0)
           monteq=min(mont(0,nonz,thchoice))
       endif
       stream=dblarr(nlon,nlatwa)
       stream(*,*)=0.
       for j=0,nlatwa-1 do begin
           nonz=where(mont(*,j,thchoice) gt 0)
           if nonz(0) ne -1 then begin
               stream(nonz,j)=crel(j) $
                         +(mont(nonz,j,thchoice)-monteq)/(2*omega*mu(j))
           endif
       endfor
       maxmont=max(stream,min=minmont)
       print,'streamfunction max, min',maxmont,minmont
       dmont=(maxmont-minmont)/(nclev-2)
       attrchoice=iv
       units='streamfunction  (m!e2!ns!e-1!n)'
       lev0=-1.00 & dlev=0.25 & llog=0 & lfill=1 & nozer=0
       lsurf=0 & lisen=1
       plotxy,idatadate,qth*pvu,longitude,latwa,alev $
         ,lsurf,lisen,thlev,attrchoice,units $
         ,lev0,dlev,nclev,llog,lfill,nozer,thchoice,projchoice
       lev0=minmont-2*dmont & dlev=dmont & llog=0 & lfill=0 & nozer=0
       lsurf=1 & lisen=1
       plotxy,idatadate,stream,longitude,latwa,alev $
         ,lsurf,lisen,thlev,attrchoice,units $
         ,lev0,dlev,nclev,llog,lfill,nozer,thchoice,projchoice
       goto,skip
   endif

   pzon=p00*total(ptop(*,0:nlatwa-1),1)/float(nlon)
   pzonzm=pzon
   qnear=0.7*maxq
   highq=where(qvert gt qnear)
   nhq=n_elements(highq)
   uave=total(utop(highq))/float(nhq)
   print,' Average u where q > Qnear = ',uave
   if lmode eq 1 and uave gt 0. and nlat gt nlatnh then begin
       print,' Winter polar vortex exists centred on max PV'
       lonc=lonc*!dpi/180.
       latc=latc*!dpi/180.
       xc=cos(lat2d)*cos(lon2d)
       yc=cos(lat2d)*sin(lon2d)
       zc=sin(lat2d)
       dphi=!dpi*0.5-latc
       xt=cos(dphi)*cos(lonc)*xc+cos(dphi)*sin(lonc)*yc-sin(dphi)*zc
       yt=-sin(lonc)*xc+cos(lonc)*yc
       zt=sin(dphi)*cos(lonc)*xc+sin(dphi)*sin(lonc)*yc+cos(dphi)*zc
       outob=where(abs(zt) gt 1.)
       if outob(0) ne -1 then begin
           zt(outob)=zt(outob)/abs(zt(outob))
       endif
       latt=asin(zt)
       xt=xt/cos(latt)
       outob=where(abs(xt) gt 1.)
       if outob(0) ne -1 then begin
           xt(outob)=xt(outob)/abs(xt(outob))
       endif
       lont=acos(xt)
       sint=yt/cos(latt)
       iwrap=where(sint lt 0.)
       lont(iwrap)=2*!dpi-lont(iwrap)
       lont=lont*180./!dpi
       latt=latt*180./!dpi
       nrlat=where(latitude lt latc*180/!dpi)
       nrlat=nrlat(0)
       if lpvplot eq 3 then begin
           plots,lonc*180/!dpi,latc*180/!dpi,psym=4,symsize=2
           plots,lont(*,nrlat),latt(*,nrlat)
           plots,lon2d(*,nrlat)*180/!dpi,lat2d(*,nrlat)*180/!dpi
       endif
;
;  Estimate zonal averages in transformed coordinates without 
;  re-gridding data which would be very slow.
;
       maxlat=90.
       minlat=latitude(1)
       inband=where(latt le maxlat and latt gt minlat)
       if inband(0) ne -1 then begin
           nband=n_elements(inband)
           pzon(0)=p00*total(ptop(inband))/float(nband)
       endif
       for j=1,nlatwa-1 do begin
           maxlat=latitude(j-1)
           minlat=latitude(j+1)
           inband=where(latt le maxlat and latt gt minlat)
           if inband(0) ne -1 then begin
               nband=n_elements(inband)
               pzon(j)=p00*total(ptop(inband))/float(nband)
           endif
       endfor
   endif else begin
;
;  When polar vortex does not exist (e.g., summer) just use 
;  zonal average p-top without transforming coordinates.
;
       print,' Not westerly stratospheric polar vortex conditions'
       print,' Use zonal average for p-top'
   endelse
   if lpvplot eq 3 then begin
       stop
       plot,latitude(0:nlatwa-1),pzon $
         ,title=' Zonal average pressure (Pa)' $
         ,xstyle=1,xrange=[0,90],xtitle='Latitude (deg)'
       plots,latitude(0:nlatwa-1),pzonzm,linestyle=1
   endif
;
;  Pre-process thetalb to remove high values that occupy hardly any 
;  surface area (i.e., mountains).
;  Work downwards in theta-levels.
;
;  Introduces problems below ground under "hot spots" because
;  effectively dropping lower boundary below ground.
;
;   nabthres=nlon
;   for m=nthlev-1,0,-1 do begin
;       thm=thlev(m)
;       above=where(thetalb gt thm)
;       nab=n_elements(above)
;       if nab lt nabthres and above(0) ne -1 then begin
;           print,'Removing surface hot spots for theta > ',thm,m,nab
;           thetalb(above)=thm
;       endif
;   endfor
;
   if lmode eq 1 or lmode eq 4 then begin
;
;  Calculation of area, mass and circulation integrals 
;  within (modified) PV contours in isentropic layers.
;  Also calculate the integrals for isentropic polar shells.
;  Choose between global or Northern hemisphere only for LB theta.
;
       lbox=boxarea(*,0:nlatend)
       thsmax=max(thetalb(*,0:nlatwa-1),min=thsmin)
       print,' Integration for area, mass and circulation'
       totmassst=0.
       tottracst=0.
       for m=0,nthlev-1 do begin
           thm=thlev(m)
;           thmmh=thlevh(m)
;           thmph=thlevh(m+1)
           thmmh=thm
           thmph=thm
;
;  Integrating over polar cap shells for layers intersecting ground.
;  Allowing for fractional layer contributions near intersection.
;
           interior=where(thetalb(*,0:nlatwa-1) lt thmmh) 
           if interior(0) gt -1 then begin
               lmass=denth(*,0:nlatwa-1,m)
               lcirc=zetamod(*,0:nlatwa-1,m)
               psarea(m)=total(lbox(interior))
               psmass(m)=total(lmass(interior)*lbox(interior))
               pscirc(m)=total(lcirc(interior)*lbox(interior))
;               psvmom2(m)=total((lcirc(interior)^2)*lbox(interior))
;               psvmom3(m)=total((lcirc(interior)^3)*lbox(interior))
;               psvmom4(m)=total((lcirc(interior)^4)*lbox(interior))
           endif           
           edge=where(thetalb(*,0:nlatwa-1) ge thmmh $
                  and thetalb(*,0:nlatwa-1) lt thmph)
           if edge(0) gt -1 then begin
               stop
               thfrac=(thmph-thetalb(*,0:nlatwa-1))/(thmph-thmmh)
               larea=lbox(edge)*thfrac(edge)
               denthm=denths(*,0:nlatwa-1)
               lmass=denthm(edge)*larea
               zetam=zetamod(*,0:nlatwa-1,m)
               lcirc=zetam(edge)*larea
               psarea(m)=psarea(m)+total(larea)
               psmass(m)=psmass(m)+total(lmass)
               pscirc(m)=pscirc(m)+total(lcirc)
           endif
;
;  Check that isentropic layer shells occupy more than a minimum area.
;  Otherwise put them below ground.
;
           if psarea(m) lt areacut*0.1 then begin
               print,' Area of surface too small and put below ground ',m,thm
               psarea(m)=0.
               psmass(m)=0.
               pscirc(m)=0.
;               psvmom2(m)=0.
;               psvmom3(m)=0.
;               psvmom4(m)=0.
               thsmin=thm+0.1
           endif
       endfor
;
;  Integrating within PV contours in each isentropic layer.
;  Set integrals for qk <= qminth to values for polar cap shells.
;
;  Loop from top of domain to bottom.
;
       for m=nthlev-1,0,-1 do begin
           thm=thlev(m)
           dth=thlevh(m+1)-thlevh(m)
           if thm ge thsmin then begin
              qmaxth(m)=max(qmod(*,0:nlatwa-1,m),min=qmin)
              qminth(m)=qmin
              print,' Level ',m,thlev(m),qmaxth(m),qminth(m)
              qsmall=where(trlevst le qmin)
              nsmall=n_elements(qsmall)
              if qsmall(0) eq -1 then begin
                 qsmall=0
              endif
              lmass=denth(*,0:nlatend-1,m)
              lcirc=zetamod(*,0:nlatend-1,m)
              klast=0
              for k=nsmall,ntrlevst-1 do begin
                 mpv=trlevst(k)
                 highq=where(qmod(*,0:nlatend-1,m) gt mpv)
                 if highq(0) gt -1 then begin
                    bsarea(k,m)=total(lbox(highq))
                    bsmass(k,m)=total(lmass(highq)*lbox(highq))
                    bscirc(k,m)=total(lcirc(highq)*lbox(highq))
                    klast=k
                 endif           
              endfor
;
; Account for a maximum in circulation which only occurs when
; negative PV is included within the domain of integration.
; The inverter cannot cope with -ve f*PV so modify the integrals.
;
              circmax=max(bscirc(*,m))
              if lnhonly eq 1 and circmax gt pscirc(m) then begin
                  highc=where(bscirc(*,m) gt pscirc(m))
                  nc=n_elements(highc)
                  topc=highc(nc-1)+1
                  if thm ge thsmax then begin
;
; For isentropic surfaces that do not intersect the ground.
; The whole region where circ > circ_NH is modified.
;
;                      print,'Modifying integrals for surface ',m,thlev(m) $
;                        ,' k = ',highc

;                      thtest=thsmax ; to use zero relative vorticity
                      thtest=thtop ; to use zero PV across region
                      if thm le thtest then begin
;
; Zero PV and absolute vorticity modification
; 
                          bsarea(highc,m)=bsarea(topc,m)
                          bsmass(highc,m)=bsmass(topc,m)
                          bscirc(highc,m)=bscirc(topc,m)
                      endif else begin
;
; Zero relative vorticity modification
;
                          emut=1-2*bsarea(topc,m)
                          dentopc=(psmass(m)-bsmass(topc,m))*2/emut
                          emuc=trlevst(highc)*lait2pv(m)*dentopc*0.5/omega
                          print,'emut = ',emut
                          print,'emuc = ',emuc
                          bsarea(highc,m)=bsarea(topc,m)+0.5*(emut-emuc)
                          bsmass(highc,m)=bsmass(topc,m) $
                                         +dentopc*0.5*(emut-emuc)
                          bscirc(highc,m)=bscirc(topc,m) $
                                         +0.5*omega*(emut^2-emuc^2)
                      endelse
                  endif else begin
;
; On isentropic surfaces that intersect the ground, assign the 
; zonal wind from the highest PV contour where circ > circ_s to 
; all PV levels below this and use it to determine circulation there.
;
;                      print,'Modifying integrals for surface ',m,thlev(m) $
;                        ,' k = ',highc
;                      emu1=1-2*bsarea(highc,m)
;                      coss1=1.-emu1*emu1
;                      u1=2.*radea*bscirc(highc,m)-omega*radea*coss1
;                      print,'emu1 = ',emu1
;                      print,'u1 = ',u1
;                      bscirc(highc,m)=(u1(nc-1)+omega*radea*coss1) $
;                        /(2*radea)
;                      emu1=1-2*psarea(m)
;                      coss1=1.-emu1*emu1
;                      pscirc(m)=(u1(nc-1)+omega*radea*coss1)/(2*radea)
                  endelse
              endif
;
;  Set the integrals for smallest PV value to the integral over the 
;  entire isentropic shell (psmass) or the integral just 
;  calculated for PV level exceeding Qmin, if this integral is larger.
;  psmass can only be smaller if some negative PV on the theta-level has 
;  contributed to the integral over the whole level.
;
              bsarea(qsmall,m)=max([psarea(m),bsarea(nsmall,m)])
              bsmass(qsmall,m)=max([psmass(m),bsmass(nsmall,m)])
              bscirc(qsmall,m)=max([pscirc(m),bscirc(nsmall,m)])
;
;  If several PV contours share the same equivalent latitude next to the
;  North Pole, then chop off the highest PV until one contour remains.
;
;  If several PV contours are so close the North Pole that area <
;  areacut, then chop off the highest PV until one contour remains.
;
              if klast gt 1 then begin
                 acut=areacut
                 bsran=bsarea(0:klast,m)
                 if m le mmed then begin
                     if thm le chopths then begin
                         acut=psarea(m)*0.5
                     endif
; NOTE: or statement below implies that PV is capped at value pvtrop
;       anywhere within polar LT where m le mmed (mu > mumod).
                     ktest=where(bsran gt 0. and (bsran lt acut $
                           or trlevst(0:klast) ge pvtrop))
                 endif else begin
                     if thm gt thref+20 then begin
                         acut=0.3*areacut
                     endif
                     ktest=where(bsran gt 0. and bsran lt acut)
                 endelse
                 kcut=ktest(0)
                 ntest=n_elements(ktest)
                 if ntest gt 1 then begin
                    kover=ktest(1:*)
                    kcut=kover(0)
                    print,' cutting off high PV spike',m,thlev(m),kcut,klast
                    qmaxth(m)=0.9*trlevst(kcut)+0.1*trlevst(kcut-1)
                    circcut=trlevst(kcut)*bsmass(kcut,m)
;
;  Correct all circulation integrals for removing highest PV values
;  assuming that the mass integrals are unchanged.
;  Also correct area integrals, either assuming that relative circulation
;  (zonal flow) does not change, or that adjusted zonal flow at Q_c is
;  zero. 
;  
;  Used for first guess equivalent latitudes.
;
                    chcirc=bscirc(kcut,m)-circcut
                    bscirc(0:kcut-1,m)=bscirc(0:kcut-1,m)-chcirc
;                        muesq=(1-2.*bsarea(0:kcut-1,m))^2
;                        muesq=chcirc*2./omega+muesq
;                        bsarea(0:kcut-1,m)=0.5*(1-sqrt(muesq))
                    mue=1-2.*bsarea(kcut,m)
                    muehat=sqrt(chcirc*2./omega+mue*mue)
                    if muehat gt 1. then begin
                       muehat=1.
                    endif
                    darea=-0.5*(muehat-mue)
                    bsarea(0:kcut-1,m)=bsarea(0:kcut-1,m)+darea
;                        outran=where(muesq ge 1.)
;                        if outran(0) ne -1 then begin
;                            print,'         muesq ge 1 ',m,thlev(m)
;                            print,outran
;                            bsarea(outran,m)=0.
;                            bsmass(outran,m)=0.
;                            bscirc(outran,m)=0.
;                            kcut=outran(0)
;                            qmaxth(m)=0.9*trlevst(kcut)+0.1*trlevst(kcut-1)
;                        endif
                    bsarea(0,m)=psarea(m)
                    bscirc(0,m)=pscirc(m)
                    bsarea(kover,m)=0.
                    bsmass(kover,m)=0.
                    bscirc(kover,m)=0.
                 endif
              endif   
;
;  Set background zonal wind at intersection with ground equal to
;  value on first interior point mu(Q_k) and adjust circulation accordingly.
;  This is used because finite difference estimate of relative
;  vorticity next to intersection with the ground is poor.
;
;              if thm lt thsmax then begin
              if thm lt -1 then begin  ; this modification is not being used!
                  acut=psarea(m)
                  highc=where(bsarea(*,m) lt acut and bsarea(*,m) gt 0.)
                  nc=n_elements(highc)
                  if nc gt 1 then begin
                      emu1=1-2*bsarea(highc,m)
                      coss1=1.-emu1*emu1
                      u1=2.*radea*bscirc(highc,m)-omega*radea*coss1
                      emus1=1-2*psarea(m)
                      cosss=1.-emus1*emus1
                      us=2.*radea*pscirc(m)-omega*radea*cosss
                      print,'us, min(u1) on theta surface ',m,us,u1(1)
                      bscirc(highc(0),m)=(u1(1)+omega*radea*coss1(0))$
                        /(2*radea)
                      pscirc(m)=(u1(1)+omega*radea*cosss)/(2*radea)
                  endif
              endif
          endif
;
;  Calculate total mass and circulation in isentropic layers.
;
          totmassst=totmassst+psmass(m)*4.*!dpi*radea*radea*dth
          tottracst=tottracst+pscirc(m)*4.*!dpi*radea*radea*dth
      endfor
      stratmass=total(pzon*dmu(0:nlatwa-1))
      pstot=total(pszav(0:nlatwa-1)*dmu(0:nlatwa-1))
      if lbound eq 1 then begin
          pt=alev(nlev-nbound)*pstot*p00
      endif else begin
          pt=p00
      endelse
      massnh=(pt-stratmass)*2*!dpi*radea*radea/ga
      print,' '
      print,'Mass from psmass / mass from psurf = ',totmassst/massnh
;
;  Calculate background zonal flow at equator from circulation
;  integrated over entire NH. 
;  Note: thlev here ordered from ground to top.
;
       above=where(thlev gt thsmax)
       ueq(above)=2.*radea*pscirc(above)-omega*radea
       print,'above = '
       print,above
       nab=above(0)
       ueq(0:nab-1)=ueq(nab)
;
;  Define orography of background state to be zonally averaged 
;  position of the lower boundary.
;
       zs0=total(zlb(*,0:nlatwa-1),1)/float(nlon)
;
;  Apply running mean smoothing to orography to avoid introducing
;  noise into the inverter solution.
;
       boxw=10
       zs0=smooth(zs0,boxw)
;
;  First guess equivalent latitudes.
;
       emu=1.-2.*bsarea
       equivlat=asin(emu)
       equivcossq=(1.-emu^2)
       equivcos=sqrt(equivcossq)
       emus=1.-2.*psarea
       cosm=sqrt(1.-emus^2)
;
; Reverse theta levels so that they run from top to bottom for
; compatibility with PV inverter.
;
       qmaxth=reverse(qmaxth) 
       qminth=reverse(qminth)  
       rlait2pv=reverse(lait2pv)  
       bsarea=reverse(bsarea,2)
       bsmass=reverse(bsmass,2)
       bscirc=reverse(bscirc,2)
       psarea=reverse(psarea)
       psmass=reverse(psmass)
       pscirc=reverse(pscirc)
       ueq=reverse(ueq)
       bsout=fpathbs+'bs_'+fstembi+strtrim(idatadate,2)
       print,'WRITING '+bsout
       openw,2,bsout
       printf,2,' DATA BASE TIME IS ',idatadate,format='(A19,I10)'
       printf,2,nlatwa,' LATITUDES ON GAUSSIAN GRID AT',format='(I4,A)'
       printf,2,latwa,format='(8E13.5)'
       printf,2,nthlev,' ISENTROPIC LEVELS AT',format='(I4,A)'
       printf,2,thlevrev,format='(10F8.2)'
       printf,2,thtop,topminth $
        ,' TOP BOUNDARY IN ISENTROPIC COORDS, min(THETA(L1))',format='(2F8.2,A)'
       printf,2,thsmax,thsmin,thoceanmax,' MAX, MIN AND EQUATORIAL VALUES OF SURFACE THETA' $
         ,format='(3E13.5,A)'
       printf,2,ntrlevst,' TRACER MIXING RATIO CONTOURS AT',format='(I4,A)'
       printf,2,trlevst,format='(8E13.5)'
       printf,2,' MAX PV ON THETA LEVELS',format='(A)'
       printf,2,qmaxth,format='(8E13.5)'
       printf,2,' MIN PV ON THETA LEVELS',format='(A)'
       printf,2,qminth,format='(8E13.5)'
       printf,2,' FACTOR TO CONVERT FROM LAIT TO ERTEL PV',format='(A)'
       printf,2,rlait2pv,format='(8E13.5)'
       printf,2,' TOTAL PVS ENCLOSED BY LOWEST VALUE TRACER CONTOUR'
       printf,2,tottracst,format='(8E13.5)' 
       printf,2,' TOTAL ATM MASS ENCLOSED BY LOWEST VALUE TRACER CONTOUR'
       printf,2,totmassst,format='(8E13.5)'
       printf,2,' AREA INTEGRALS IN PV-THETA COORDINATES' 
       printf,2,bsarea,format='(8E13.5)'
       printf,2,' MASS INTEGRALS IN PV-THETA COORDINATES' 
       printf,2,bsmass,format='(8E13.5)'
       printf,2,' CIRCULATION INTEGRALS IN PV-THETA COORDINATES' 
       printf,2,bscirc,format='(8E13.5)'
       printf,2,' AREA INTEGRAL OVER POLAR SHELLS THETA COORDINATES' 
       printf,2,psarea,format='(8E13.5)'
       printf,2,' MASS INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES' 
       printf,2,psmass,format='(8E13.5)'
       printf,2,' CIRCULATION INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES'
       printf,2,pscirc,format='(8E13.5)'
       printf,2,' BACKGROUND PRESSURE ON TOP BOUNDARY'
       printf,2,pzon,format='(8E13.5)'
       printf,2,' BACKGROUND SURFACE GEOPOTENTIAL'
       printf,2,zs0,format='(8E13.5)'
       printf,2,' BACKGROUND u*cos(phi) AT EQUATOR ON THETA LEVELS'
       printf,2,ueq,format='(8E13.5)'
       printf,2,' max Lait PV above stratopause and average u (q > 0.7 Q_max)'
       printf,2,maxq,uave,format='(8E13.5)'
       close,2

       lmoments=0
       if lmoments eq 1 then begin
           psvmom2=reverse(psvmom2)
           psvmom3=reverse(psvmom3)
           psvmom4=reverse(psvmom4)
           vsout='vortmoments_'+fstembi+strtrim(idatadate,2)
           print,'WRITING '+vsout
           openw,2,vsout
           printf,2,' DATA BASE TIME IS ',idatadate,format='(A19,I10)'
           printf,2,nlatwa,' LATITUDES ON GAUSSIAN GRID AT',format='(I4,A)'
           printf,2,latwa,format='(8E13.5)'
           printf,2,nthlev,' ISENTROPIC LEVELS AT',format='(I4,A)'
           printf,2,thlevrev,format='(10F8.2)'
           printf,2,' AREA INTEGRAL OVER POLAR SHELLS THETA COORDINATES' 
           printf,2,psarea,format='(8E13.5)'
           printf,2,' MASS INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES' 
           printf,2,psmass,format='(8E13.5)'
           printf,2,' CIRCULATION INTEGRALS OVER POLAR SHELLS IN THETA COORDINATES'
           printf,2,pscirc,format='(8E13.5)'
           printf,2,' SECOND MOMENT OF ABSOLUTE VORTICITY IN THETA COORDINATES'
           printf,2,psvmom2,format='(8E13.5)'
           printf,2,' THIRD MOMENT OF ABSOLUTE VORTICITY IN THETA COORDINATES'
           printf,2,psvmom3,format='(8E13.5)'
           printf,2,' FOURTH MOMENT OF ABSOLUTE VORTICITY IN THETA COORDINATES'
           printf,2,psvmom4,format='(8E13.5)'
           close,2
       endif
   endif else begin
;
;  Read integral diagnostics output by mode 1.
;
       if lbsevol eq 1 then begin
           bsname=fpathbs+'bs_'+fstembi+strtrim(idatadate,2)
       endif else begin
           bsname=fpathbs+'bs_'+fstembi+'2001010000'
       endelse
       print,'READING ',bsname
       openr,2,bsname
       readf,2,scrap
       readf,2,nlatin,scrap,format='(I4,A)'
       if nlatin ne nlatwa then begin
           print,'nlatin ne nlatwa ',nlatin,nlatwa
           stop
       endif
       latin=dblarr(nlatin)
       readf,2,latin,format='(8E13.5)'
       readf,2,nthlevin,scrap,format='(I4,A)'
       if nthlevin ne nthlev then begin
           print,'nthlevin ne nthlev ',nthlevin,nthlev
           stop
       endif
       thlevin=thlev
;       readf,2,thlevin,format='(10F8.3)'
;       readf,2,thtop,scrap,format='(F8.3,A)'
       readf,2,thlevin,format='(10F8.2)'
       readf,2,thtop,scrap,format='(F8.2,A)'
       print,thtop,scrap
       readf,2,thsmax,thsmin,scrap,format='(2E13.5,A)'
       print,thsmax,thsmin,scrap
       readf,2,ntrlevin,scrap,format='(I4,A)'
       if ntrlevin ne ntrlevst then begin
           print,'ntrlevin ne ntrlevst ',ntrlevin,ntrlevst
           stop
       endif
       trlevin=trlevst
       readf,2,trlevin,format='(8E13.5)'
       readf,2,scrap
       readf,2,qmaxth,format='(8E13.5)'
       readf,2,scrap
       readf,2,qminth,format='(8E13.5)'
       if llait eq 1 then begin
           readf,2,scrap
           readf,2,lait2pv,format='(8E13.5)'
       endif
       readf,2,scrap
       readf,2,tottracst,format='(8E13.5)' 
       readf,2,scrap
       readf,2,totmassst,format='(8E13.5)'
       readf,2,scrap
       readf,2,bsarea,format='(8E13.5)'
       readf,2,scrap
       readf,2,bsmass,format='(8E13.5)'
       readf,2,scrap
       readf,2,bscirc,format='(8E13.5)'
       readf,2,scrap
       readf,2,psarea,format='(8E13.5)'
       readf,2,scrap
       readf,2,psmass,format='(8E13.5)'
       readf,2,scrap
       readf,2,pscirc,format='(8E13.5)'
       readf,2,scrap
       readf,2,pzon,format='(8E13.5)'
       readf,2,scrap
       readf,2,zs0,format='(8E13.5)'
       readf,2,scrap
       readf,2,ueq,format='(8E13.5)'
       close,2
       lait2pv=reverse(lait2pv)
   endelse
;
; Re-order arrays so that theta runs from bottom to top.
;
   qmaxth=reverse(qmaxth)
   qminth=reverse(qminth)
   bsarea=reverse(bsarea,2)
   bsmass=reverse(bsmass,2)
   bscirc=reverse(bscirc,2)
   psarea=reverse(psarea)
   psmass=reverse(psmass)
   pscirc=reverse(pscirc)
   ueq=reverse(ueq)
;
;  Spawn script to run PV inverter using bsout file (written above) as input.
;
   if lmode eq 4 then begin
       openw,2,fpathjob+'current.dat'
       printf,2,idatadate,format='(i10)'
       close,2
       scrnam=fpathjob+'zminv.bl'
       outnam=fpathjob+'out'+strtrim(idatadate,2)
       command=scrnam+' > '+outnam
       print,'SPAWNING PV INVERTER SCRIPT '+command
       spawn, /sh, command
       print,'RETURNED from '+command
   endif
   if lmode gt 1 then begin
       ;emus(*)=0.
       ;emu(*,*)=0.
       ; *****************************************************************
       ; Load background state data
       ;bs_file = fpathpv+'bs_data_'+string(idatadate, format='(I010)')+'.nc'
       ;id = NCDF_OPEN(bs_file)
      
       ; Read equivalent latitudes output by outer iteration of PV inverter.
       ; equivalent latitudes (PV-theta grid)
       ;emu_id = NCDF_VARID(id, 'emu')
       ;NCDF_VARGET, id, emu_id, emu
       ; surface equivalent latitudes (theta grid)
       ;emus_id = NCDF_VARID(id, 'surf_emu')
       ;NCDF_VARGET, id, emus_id, emus
       ; Re-order arrays so that theta runs from bottom to top.
       ;emu=reverse(emu,2)
       ;emus=reverse(emus)
       ;
       ; Ensure monotonicity of equivalent latitudes wrt PV. 
       ;for l=0,nthlev-1 do begin
       ;    minemu=min(emu(*,l),minel)
       ;    if minel gt 0 then begin
       ;        print,'Forcing monotonic emu on theta = ',thlev(l),minel
       ;        emu(0:minel-1,l)=emus(l)
       ;    endif
       ;endfor
       ; *****************************************************************
      
       equivlat=asin(emu)
       equivcossq=(1.-emu^2)
       equivcos=sqrt(equivcossq)
       cosm=sqrt(1.-emus^2)
       dfac=4.*!dpi*radea*radea/(2.*!dpi*radea)
;
;  Note that pvinv contains Ertel PV and has been transformed to
;  pv0 which contains modified (Lait) PV and is ordered from NP to eq.
;
;  Calculate the background zonal flow at surface lats mu_es(theta_m).
;
       bsucosm=2.*radea*pscirc-radea*omega*(1.-emus*emus)
       notnp=where(emus lt 1.)
       bsucosm(notnp)=bsucosm(notnp)/cosm(notnp)
;
; Interpolations needed to map background state variables from
; (Q_k, theta_m) to (mu_j, theta_m) using equivalent latitudes.
;
       wa3d(*,*,*)=0.
       psik=dblarr(ntrlevst,nthlev)
       psik(*,*)=0.
       ecask=psik
       psi=dblarr(ntrlevst)
;
; First, linearly interpolate surface u from grid to equivalent latitudes
; mu_es(theta_m) - to be used as a boundary condition.
;
       usm=bsucosm
       useq=0.
       intj2ths,thlev,ts0,us0,useq,thsmin,thsmax,emus,usm
;
; Now loop over theta-levels
;
       for l=0,nthlev-1 do begin
           thl=thlev(l)
           thfac=lait2pv(l)
           emumin=emu(0,l)
           emu1=emu(*,l)
           karr=where(emu1 lt 1. and (emu1-emumin) gt tiny)
           nk=n_elements(karr)
           if nk le 1 then begin
               n=-1
           endif else begin
               kmin=karr(1)
               kmax=karr(nk-1)
               if emus(l) gt 0. then begin
                   if kmin eq 0 then begin
;                       print,'kmin=1 not allowed for isentropic layer intersecting ground, since q_surf > 0 required'
                       kmin=1
                   endif
               endif
;               print,'kmin = ',kmin,'  kmax = ',kmax
               n=kmax-kmin
           endelse
;           print,n,' PV values on isentropic surface ',thl
           if n ge 0 then begin
               qedge=qminth(l)
               if qedge ge 0 then begin
                   if emus(l) eq 0 then begin
                       qedge=0.
                   endif else begin
                       qedge=trlevst(kmin-1)
                   endelse
               endif
;
;  Calculate the background state PV gradient from
;  data in (mu_j, theta_m) coordinates and its reciprocal.
;
               for j=1,nlatwa-2 do begin
                  qbaryj(j,l)=(pv0(j-1,l)-pv0(j+1,l))/(mu(j-1)-mu(j+1))
               endfor
               qbaryj(*,l)=qbaryj(*,l)/radea
               nonzero=where(qbaryj(*,l) gt 0)
               grlatj(nonzero,l)=1./qbaryj(nonzero,l)
               qbaryj(*,l)=1.e11*(1.-mu^2)*sigma0(*,l)*qbaryj(*,l)
;
;  Overwrite mass and circulation integral arrays
;  using PV derived from first guess latitudes, as done by BSMODIFY.
;
;  JM 31/7/25: Overwriting bsmass and bscirc is not needed for OT MLM state.
;
               emu1=emu(*,l)
               emusurf=emu1(0)
               pvm=pv0(*,l)
               msurf=psmass(l)
               ;intj2pv,l,thl,trlevst,pvm,mlat,qmaxth,qminth $
               ;       ,emu1,emusurf,bsmass,msurf,0.
               clat=0.5*uth0(*,l)*cosj/radea+0.5*omega*cosj*cosj
               csurf=pscirc(l)
               ;intj2pv,l,thl,trlevst,pvm,clat,qmaxth,qminth $
               ;       ,emu1,emusurf,bscirc,csurf,0.
               pcask(*,l)=2.*radea*(-bscirc(*,l) $
                                    +trlevst*thfac*(bsmass(*,l)-psmass(l)))
;
;  Interpolate zonal wind and density from inverter grid to 
;  final equivalent latitudes, using final PV from inverter.
;
               emu1=emu(*,l)
               emusurf=emu1(0)
               pvm=pv0(*,l)
               uinput=uth0(*,l)
               intj2pv,l,thl,trlevst,pvm,uinput,qmaxth,qminth $
                      ,emu1,emusurf,uthk,bsucosm(l),0.
               sinput=sigma0(*,l)
               snonz=where(sinput gt 0.)
               snel=n_elements(snonz)
               snp=sigma0(snonz(0),l)
               ssurf=sigma0(snonz(snel-1),l)
               intj2pv,l,thl,trlevst,pvm,sinput,qmaxth,qminth $
                      ,emu1,emusurf,sigmak,ssurf,snp
;
;  Integrate zonal wind to find mass streamfunction as a function of 
;  PV and theta.
;  (note uth0 and uthk are not weighted by cos(lat))
;  
;  Also, find Casimir defined by integral of background streamfunction
;  with respect to PV.
;
               psi(*)=0.
               ecas=psi
               latmin=asin(emumin)
               jels=where(rlatinv*!pi/180. gt latmin)
               jels=jels(0)
               if emus(l) eq 0. then begin
;
;  Use equatorial zonal wind from inverter.
;
                   qedge=0.
                   fedge=radea*sigmainv(0,l)*ebu0(l)
                   psiedge=0.
               endif else begin
;
;  Assign boundary streamfunction such that boundary terms are second order.
;
                   qedge=trlevst(kmin-1)*thfac
                   sedge=sigmak(kmin-1,l)
                   fedge=radea*sedge*usm(l)
;                   psiedge=ga*zsinv(jels)/qedge ; appropriate BC
                   psiedge=0. ; zero assumed for linearisation of H_e in paper
               endelse
               psi(*)=psiedge
               x0=rlatinv*!pi/180.
               f0=radea*sigmainv(*,l)*uthinv(*,l)
               for k=kmin,kmax do begin
                   latk=asin(emu(k,l))
                   muel=where(x0 lt latk and x0 gt latmin)
                   if muel(0) ne -1 then begin
                       x=[latmin,x0(muel),latk]
                       f=-[fedge,f0(muel),radea*sigmak(k,l)*uthk(k,l)]
                       result=int_tabulated(x,f,/double)
                       psi(k)=psiedge+result
                   endif
               endfor
               if kmax+1 lt ntrlevst then begin
                   latnp=!dpi/2.
                   muel=where(x0 gt latmin)
                   if muel(0) ne -1 then begin
                       x=[latmin,x0(muel),latnp]
                       f=-[fedge,f0(muel),0.]
                       result=int_tabulated(x,f,/double)
                       psi(kmax+1:*)=psiedge+result
                   endif
               endif
               for k=1,ntrlevst-1 do begin
                   qx=[trlevst(0:k)]*thfac
                   psix=[psi(0:k)]
                   result=trap_uneven(qx,psix)
                   ecas(k)=result
               endfor
               psik(*,l)=psi
               ecask(*,l)=ecas-ecas(kmin-1)
;               ecask(*,l)=ecas ; edit 13/1/16 (better cancellation in H_e)
           endif
       endfor
;
; **************************************************************************
; Calculate finite amplitude wave activity (pseudomomentum and pseudoenergy)
; Limit calculation to range of theta levels within input data (0:nthlim-1).
;
       for l=0,nthlim-1 do begin
           thl=thlev(l)
           thfac=lait2pv(l)
           psi0b=0.
           m0b=0.
           e0s=0.
           ke0s=0.
           emumin=emu(0,l)
           karr=where(emu(*,l) lt 1. and (emu(*,l)-emumin) gt tiny)
           nk=n_elements(karr)
;
; Check that thlev(l) is somewhere above ground, because occasionally
; emus(l) < 1 even though the surface is lower than the boundary theta. 
;
           smallth=where(thetalb(*,0:nlatwa-1) lt thl)
           if smallth(0) eq -1 and nk gt 1 then begin
               print,nk, $
               ' points with emu < 1 when surface lower than thetalb ',thl
               emu(*,l)=1.
               emus(l)=1.
               nk=0
           endif
           if nk le 1 then begin
               n=-1
           endif else begin
               kmin=karr(1)
               kmax=karr(nk-1)
               if emus(l) gt 0. then begin
                   if kmin eq 0 then begin
;                       print,'kmin=1 not allowed for isentropic layer intersecting ground, since q_surf > 0 required'
                       kmin=1
                   endif
               endif
;               print,'kmin = ',kmin,'  kmax = ',kmax
               n=kmax-kmin
           endelse
           print,n,' PV values on isentropic surface ',thl
           if n ge 0 then begin
               jelmax=nlatwa
               jelb=nlatwa
               if emus(l) gt 0. then begin
;
; Calculate boundary, then interior contributions to wave activity density.
; First specify bounds of integration if theta surface intersects the
; ground.
; Remember latitudes ordered from NP southwards.
;
; jels = index for latitude immediately south of bg state LB (emus)
; jelmin = index for latitude of most southern excursion of theta_l
; jelmax = index for latitude of most northern excursion of theta_l
; Note latitude ordered from NP towards equator so we expect
; jelmax < jels < jelmin
;
                   jels=where(mu(0:nlatwa-1) le emus(l)) & jels=jels(0)
                   mumin=min(muarr(smallth))
                   jelmin=where(mu(0:nlatwa-1) eq mumin) & jelmin=jelmin(0)
                   largeth=where(thetalb(*,0:nlatwa-1) ge thl)
                   if largeth(0) gt -1 then begin
                       mumax=max(muarr(largeth))
                       jelmax=where(mu(0:nlatwa-1) eq mumax) & jelmax=jelmax(0)
                   endif else begin
;
; Theta surface does not intersect ground in 3D state.
;
                       mumax=emus(l)
                       jelmax=jels & jelmin=jels
                   endelse
                   print,l,thl,jels,emus(l),jelmax,jelmin
                   if jelmax gt jelmin then begin
;
; Prevents error where max(lat) < min(lat) for wavy state.
;
                       jelmax=jelmin
                   endif
                   if jels gt jelmin then begin
;
; Reduce basic state domain so that it does not lie completely 
; outside wavy state domain. 
;
                       jels=jelmin
                   endif
                   if jels le jelmax then begin
                       if jels eq -1 then begin
;
; Treat as if surface is entirely above ground.
;
                           jelmax=nlatwa
                           jelmin=-1
                       endif else begin
;
; Reduce wavy state domain so that it does not lie completely 
; outside basic state domain. 
;
                           jelmax=jels
                           jelmin=jels
                       endelse
                   endif
                   jelb=jelmax-1
;
; Find surface energy as boundary condition on
; integration to find energy Casimir.
;
                   if jels gt 0 then begin
                       j=jels-1
                       ke0s=0.5*us0(j)*us0(j)
                       e0s=ke0s+cp*thl*((ps0(j)/p00)^kappa)
;                       ets=ke0s+cp*thl*((plbav(j))^kappa)
                       ets=e0s
                   endif 
;
; Find densities in domain between mumax and emus(l)
; Within 2D state domain, but outside D-bar.
;                  
                   if jelmax lt jels then begin
                       for j=jels-1,jelmax,-1 do begin
                           mpvt=pv0(j,l) ; Lait (modified) PV
                           q0=mpvt*thfac ; converted to Ertel PV
                           p0=pr0(j,l)
                           cossq=cosj(j)*cosj(j)
                           planmom=omega*radea*cossq
                           den0i=sigma0(j,l)
                           u0i=uth0(j,l)
                           den0=den0i
                           u0=u0i

                           if den0 eq 0. or temp0(j,l) eq 0. then begin
                               print,'Zero BS density or T at j,l = ',j,l
                               stop
                           endif
                           e0=0.5*u0*u0+cp*temp0(j,l)
                           d2hdpdth=cp*((p0/p00)^kappa) $
                             *(kappa/p0-kappa*(kappa-1)*thl*ga*den0i/(p0*p0))
                           karr=where(trlevst gt mpvt)
                           k=karr(0)
                           if k gt 0 then begin
                               if emu(k,l) lt 1. then begin
                                   rj=(mpvt-trlevst(k-1)) $
                                     /(trlevst(k)-trlevst(k-1))
                                   m0=(1.-rj)*bsmass(k-1,l)+rj*bsmass(k,l)
                                   c0=(1.-rj)*bscirc(k-1,l)+rj*bscirc(k,l)
                                   psi0=(1.-rj)*psik(k-1,l)+rj*psik(k,l)
                                   ecas0=(1.-rj)*ecask(k-1,l)+rj*ecask(k,l)
                               endif else begin
;
;               Close to North Pole use B.C. A=M=U=0 at NP.
;
;                                   print,'In D0: conditional clause NP'
                                   rj=(mpvt-trlevst(kmax)) $
                                     /(qmaxth(l)-trlevst(kmax))
                                   rj=min([rj,1.])
                                   omrj=1.-rj
                                   emu0=omrj*emu(kmax,l)+rj*1.
                                   m0=omrj*bsmass(kmax,l)
                                   ucosj=omrj*uthk(kmax,l)*equivcos(kmax,l)
                                   c0=ucosj/(2.*radea)+0.5*omega*(1.-emu0*emu0)
                                   psi0=psik(kmax,l)
                                   ecas0=ecask(kmax,l) $
                                        +psi0*(mpvt-trlevst(kmax))*thfac
                               endelse
                           endif else begin
;                               print,'In D0: conditional clause k le 0'
                               if k eq 0 then begin
;
;               At equator.
;                                       
                                   m0=psmass(l)
                                   c0=pscirc(l)
                                   psi0=psik(0,l)
                                   ecas0=ecask(0,l)
                               endif else begin
;
;               Beyond range of PV values.
;
                                   m0=bsmass(ntrlevst-1,l)
                                   c0=bscirc(ntrlevst-1,l)
                                   psi0=psik(ntrlevst-1,l)
                                   ecas0=ecask(ntrlevst-1,l)
                               endelse
                           endelse
                           ecas0=ecas0-e0s
                           walinav=0.
                           wad1=0.
                           wad2=0.
                           wag=0.
                           wae2=0.
                           ped1=0.
                           ped2=0.
                           peke=0.
                           pep=0.
                           pes=0.
                           peg=0.
                           pee2=0.
                           fphi=0.
                           fth=0.
                           for i=0,nlon-1 do begin
                               st=denth(i,j,l)
                               if st gt 0. then begin
;
; Within fluid domain defined by 3D distribution too (i.e., intersection).
;
                                   ;rwt=0.5*(st+den0) ; density weight
                                   rwt=st
                                   mpvt=qth(i,j,l)  ; Lait PV
                                   qt=mpvt*thfac    ; Converted to Ertel PV
                                   ue=uthf(i,j,l)/cosj(j)-u0
                                   ve=vthf(i,j,l)/cosj(j)
                                   pt=p00*(tth(i,j,l)/thl)^(1./kappa)
                                   pe=pt-p0
                                   karr=where(trlevst gt mpvt)
                                   k=karr(0)
                                   if k gt 0 then begin
                                       if emu(k,l) lt 1. then begin
                                           rj=(mpvt-trlevst(k-1)) $
                                             /(trlevst(k)-trlevst(k-1))
                                           mt=(1.-rj)*bsmass(k-1,l) $
                                             +rj*bsmass(k,l)
                                           ct=(1.-rj)*bscirc(k-1,l) $
                                             +rj*bscirc(k,l)
                                           cast=(1.-rj)*pcask(k-1,l) $
                                               +rj*pcask(k,l)
                                           ecast=(1.-rj)*ecask(k-1,l) $
                                                +rj*ecask(k,l)
                                       endif else begin
;
;               Close to North Pole use B.C. A=M=U=0 at NP.
;
;                                          print,'In D0: conditional clause NP'
                                           rj=(mpvt-trlevst(kmax)) $
                                             /(qmaxth(l)-trlevst(kmax))
                                           rj=min([rj,1.])
                                           omrj=1.-rj
                                           emut=omrj*emu(kmax,l)+rj*1.
                                           mt=omrj*bsmass(kmax,l)
                                           ucosj=omrj*uthk(kmax,l) $
                                                     *equivcos(kmax,l)
                                           ct=ucosj/(2.*radea) $
                                            +0.5*omega*(1.-emut*emut)
                                           cast=omrj*pcask(kmax,l) $
                                      +rj*qmaxth(l)*thfac*2.*radea*(-psmass(l))
                                           ecast=ecask(kmax,l) $
                                      +psik(kmax,l)*(mpvt-trlevst(kmax))*thfac
                                       endelse
                                   endif else begin
;                                       print,'In D0: conditional clause k le 0'
                                       if k eq 0 then begin 
                                           mt=psmass(l)
                                           ct=pscirc(l)
                                           cast=-2.*radea*pscirc(l)
                                           ecast=ecask(0,l)
                                       endif else begin
                                           mt=bsmass(ntrlevst-1,l)
                                           ct=bscirc(ntrlevst-1,l)
                                           cast=pcask(ntrlevst-1,l)
                                           ecast=ecask(ntrlevst-1,l)
                                       endelse
                                   endelse
                                   ecast=ecast-e0s
                                   walinav=walinav+(qt-q0)^2
                                   watmp=c2fix*rwt*dfac*(-qt*(mt-m0)+ct-c0)
                                   wad1=wad1+watmp
                                   wad2=wad2-c2fix*den0*ue*cosj(j) $
                                     -c2fix*dfac*(m0-psmass(l))*(st*qt-den0*q0)
                                   wagtmp=-(st-den0)*ue*cosj(j)
                                   wag=wag+wagtmp
                                   h2=thl*cp*((pt/p00)^kappa) $
                                     -thl*cp*((p0/p00)^kappa) $
                                            *(1-kappa+kappa*pt/p0)
                                   c2=ecast-ecas0-psi0*(qt-q0)
                                   peke=peke+rwt*0.5*(ue*ue+ve*ve)
                                   pes=pes+c2fix*rwt*c2
                                   peg=peg+(st-den0)*u0*ue
                                   pep=pep+rwt*h2+0.5*pe*pe*d2hdpdth/ga
                                   ped1=ped1+c2fix*den0*u0*ue $
                                            +c2fix*psi0*(st*qt-den0*q0)
;                                   wad1=wad1-st*(uthf(i,j,l)+planmom+cast)
;                                   wad2=wad2+den0*q0*dfac*(m0-psmass(l))
;                                   ped1=ped1+st*(u0*ue+0.5*(ue*ue+ve*ve) $
;                                                +c2)+den0*h2 $
;                                            +0.5*pe*pe*d2hdpdth/ga $
;                                            +psi0*(st*qt-den0*q0)
                                   fphi=fphi+ve*(watmp+wagtmp) $
                                            -den0*ue*ve*cosj(j)
                                   fth=fth+pe*dmdlam(i,j,l)/(radea*ga)
                               endif else begin
;
; Within fluid domain defined by 2D background state, but not 3D distribution.
;
                                   wae2=wae2+den0*q0*dfac*(m0-psmass(l))
                                   pee2=pee2-den0*(e0+ecas0)
                               endelse
                           endfor
                           if ldterms eq 1 then begin
                               wad(j,l)=(wad1+wad2+wag)*dlon
                               ped(j,l)=(peke+pep+pes+ped1+peg)*dlon
;                               wad(j,l)=(wad1+wad2)*dlon
;                               ped(j,l)=(peke+pep+pes+ped1)*dlon
                           endif else begin
                               waden(j,l)=wad1*dlon
                               wagrav(j,l)=wag*dlon
                               wad(j,l)=wad2*dlon
                               peden(j,l)=peke*dlon
                               peape(j,l)=pep*dlon
                               pew(j,l)=pes*dlon
                               pegrav(j,l)=peg*dlon
                               ped(j,l)=ped1*dlon
                           endelse
                           wae(j,l)=wae2*dlon
                           waetmp(j,l)=wae2*dlon
                           pee(j,l)=pee2*dlon
                           peetmp(j,l)=pee2*dlon
                           fluxphi(j,l)=fphi*dlon
                           fluxth(j,l)=fth*dlon
                       endfor
                       jelb=jelmax
                       psi0b=psi0
                       m0b=m0
                   endif
;
; Find densities in domain between emus(l) and mumin
;                  
;                   if jels+1 lt jelmin then begin   ; introduced for ERA-I
;                       for j=jelmin,jels+1,-1 do begin
                   if jels lt jelmin then begin
                       for j=jelmin,jels,-1 do begin
                           cossq=cosj(j)*cosj(j)
                           planmom=omega*radea*cossq
                           ps0j=ps0(j)
                           wae1=0.
                           pee1=0.
                           for i=0,nlon-1 do begin
                               st=denth(i,j,l)
                               if st gt 0. then begin
;
; Within fluid domain defined by 3D distribution, but not 2D state.
;
                                   mpvt=qth(i,j,l)
                                   qt=mpvt*thfac
                                   pt=p00*(tth(i,j,l)/thl)^(1./kappa)
                                   et=0.5*(uthf(i,j,l)^2+vthf(i,j,l)^2)/cossq $
                                      +cp*tth(i,j,l)
;                                   ets=ke0s $
;                                      +cp*thetalb(i,j)*(plb(i,j)^kappa)
                                   karr=where(trlevst gt mpvt)
                                   k=karr(0)
                                   if k gt 0 then begin
                                       if emu(k,l) lt 1. then begin
                                           rj=(mpvt-trlevst(k-1)) $
                                             /(trlevst(k)-trlevst(k-1))
                                           cast=(1.-rj)*pcask(k-1,l) $
                                               +rj*pcask(k,l)
                                           ecast=(1.-rj)*ecask(k-1,l) $
                                                +rj*ecask(k,l)
                                       endif else begin
;
;               Close to North Pole use B.C. A=M=U=0 at NP.
;
;                                           print,'Out D0: conditional clause NP ',i,j,k,kmax
                                           rj=(mpvt-trlevst(kmax)) $
                                             /(qmaxth(l)-trlevst(kmax))
                                           rj=min([rj,1.])
                                           omrj=1.-rj
                                           cast=omrj*pcask(kmax,l) $
                                      +rj*qmaxth(l)*thfac*2.*radea*(-psmass(l))
                                           ecast=ecask(kmax,l) $
                                      +psik(kmax,l)*(mpvt-trlevst(kmax))*thfac
                                       endelse
                                   endif else begin
;                                       print,'Out D0: conditional clause k le 0'
                                       if k eq 0 then begin
                                           cast=-2.*radea*pscirc(l)
                                           ecast=ecask(0,l)
                                       endif else begin
;                                          Only if mpvt > trlevst(*)
                                           cast=pcask(ntrlevst-1,l)
                                           ecast=ecask(ntrlevst-1,l)
                                       endelse
                                   endelse
                                   ecast=ecast-ets
                                   wae1=wae1-st*(uthf(i,j,l)+planmom+cast)
                                   pee1=pee1+st*(et+ecast)
                               endif
                           endfor
                           wae(j,l)=wae1*dlon
                           pee(j,l)=pee1*dlon
                       endfor                   
                   endif
               endif
;
; Calculate the interior wave activity density (within D-bar).
;
               for j=jelmax-1,0,-1 do begin
                   mpvt=pv0(j,l)          ; Lait (modified) PV
                   q0=mpvt*thfac          ; converted to Ertel PV
                   dmudq0=grlatj(j,l)
                   cossq=cosj(j)*cosj(j)
                   p0=pr0(j,l)                   
                   den0i=sigma0(j,l)
                   u0i=uth0(j,l)          ; zonal wind (unweighted by cosj)
                   den0=den0i
                   u0=u0i

                   if p0 eq 0 or den0 eq 0 then begin
                       print,'Pressure/density from inverter are zero here.'
                       print,stop
                   endif
;
;                  This c0 should not be used because over-written below.
;                  Also, this formula gives the circulation of the 2-D
;                  state at lat index j, which need not match the 
;                  circulation integral at equivalent latitude mu(Q)
;                  which is estimated below by interpolation between mu(Q_k).
;
                   c0=0.5*u0*cosj(j)/radea+0.5*omega*(1.-mu(j)*mu(j))
                   e0=0.5*u0*u0+cp*temp0(j,l)
                   d2hdpdth=cp*((p0/p00)^kappa) $
                    *(kappa/p0-kappa*(kappa-1)*thl*ga*den0i/(p0*p0))
                   karr=where(trlevst gt mpvt)
                   k=karr(0)
                   if k gt 0 then begin
                       if emu(k,l) lt 1. then begin
                           rj=(mpvt-trlevst(k-1))/(trlevst(k)-trlevst(k-1))
                           m0=(1.-rj)*bsmass(k-1,l)+rj*bsmass(k,l)
                           c0=(1.-rj)*bscirc(k-1,l)+rj*bscirc(k,l)
                           psi0=(1.-rj)*psik(k-1,l)+rj*psik(k,l)
                           ecas0=(1.-rj)*ecask(k-1,l)+rj*ecask(k,l)
                       endif else begin
;
;               Close to North Pole use B.C. A=M=U=0 at NP.
;
;                           print,'In D0: conditional clause NP'
                           rj=(mpvt-trlevst(kmax)) $
                             /(qmaxth(l)-trlevst(kmax))
                           rj=min([rj,1.])
                           omrj=1.-rj
                           emu0=omrj*emu(kmax,l)+rj*1.
                           m0=omrj*bsmass(kmax,l)
                           ucosj=omrj*uthk(kmax,l)*equivcos(kmax,l)
                           c0=ucosj/(2.*radea)+0.5*omega*(1.-emu0*emu0)
                           psi0=psik(kmax,l)
                           ecas0=ecask(kmax,l)+psi0*(mpvt-trlevst(kmax))*thfac
                       endelse
                   endif else begin
                       print,'In D0: conditional clause k le 0'
                       if k eq 0 then begin
;
;               At equator.
;                                       
                           m0=psmass(l)
                           c0=pscirc(l)
                           psi0=psik(0,l)
                           ecas0=ecask(0,l)
                       endif else begin
;
;               Beyond range of PV values.
;
                           m0=bsmass(ntrlevst-1,l)
                           c0=bscirc(ntrlevst-1,l)
                           psi0=psik(ntrlevst-1,l)
                           ecas0=ecask(ntrlevst-1,l)
                       endelse
                   endelse
;
; FFT the perturbations in zonal direction to obtain Re and Im parts
; of Fourier coefficients for zonal wavenumber mselect
;
; Make sure that zonal average of filtered wavy state that has been
; obtained in isentropic coords equals MLM background state
;
                   if mselect gt 0 then begin
                       iabove=where(tth(*,j,l) gt 0.)
                       nabove=n_elements(iabove)
                       if nabove eq nlon then begin
                           upert=uthf(*,j,l)/cosj(j)
;                           upert=uthf(*,j,l)/cosj(j)-u0
                           ufour=FFT(upert,-1)
                           uifour(j,l)=ufour(mselect)
                           pmonly(mselect)=ufour(mselect)
                           pfilt=FFT(pmonly,1)
                           uav=total(upert)/nabove
                           uthf(*,j,l)=(2*real_part(pfilt)+uav)*cosj(j)
;                           uthf(*,j,l)=(2*real_part(pfilt)+u0)*cosj(j)
                           
                           vpert=vthf(*,j,l)/cosj(j)
                           vfour=FFT(vpert,-1)
                           vifour(j,l)=vfour(mselect)
                           pmonly(mselect)=vfour(mselect)
                           pfilt=FFT(pmonly,1)
                           vav=total(vpert)/nabove
                           vthf(*,j,l)=(2*real_part(pfilt)+vav)*cosj(j)
;                           vthf(*,j,l)=(2*real_part(pfilt))*cosj(j)
                           
                           ppert=p00*(tth(*,j,l)/thl)^(1./kappa)
;                           ppert=p00*(tth(*,j,l)/thl)^(1./kappa)-p0
                           pfour=FFT(ppert,-1)
                           pifour(j,l)=pfour(mselect)
                           pmonly(mselect)=pfour(mselect)
                           pfilt=FFT(pmonly,1)
                           pav=total(ppert)/nabove
                           tth(*,j,l)=thl*(((2*real_part(pfilt)+pav)/p00)^kappa)
;                           tth(*,j,l)=thl*(((2*real_part(pfilt)+p0)/p00)^kappa)
                           
                           qpert=qth(*,j,l)*thfac
;                           qpert=qth(*,j,l)*thfac-q0
                           qfour=FFT(qpert,-1)
                           qifour(j,l)=qfour(mselect)
                           pmonly(mselect)=qfour(mselect)
                           pfilt=FFT(pmonly,1)
                           qav=total(qpert)/nabove
                           qth(*,j,l)=(2*real_part(pfilt)+qav)/thfac
;                           qth(*,j,l)=(2*real_part(pfilt)+q0)/thfac
                           
                           rpert=denth(*,j,l)
;                           rpert=denth(*,j,l)-den0
                           rfour=FFT(rpert,-1)
                           pmonly(mselect)=rfour(mselect)
                           pfilt=FFT(pmonly,1)
                           rav=total(rpert)/nabove
                           denth(*,j,l)=(2*real_part(pfilt)+rav)
;                           denth(*,j,l)=(2*real_part(pfilt)+den0)
                       endif else begin
                           print,' UNEXPECTED points below ground within D-bar'
                           print,j,l,latitude(j),thlev(l)
                           denth(*,j,l)=0.
                       endelse
                   endif
                   wazonav=0.
                   walinav=0.
                   wad1=0.
                   wag=0.
                   peke=0.
                   pes=0.
                   peg=0.
                   pep=0.
                   ped1=0.
                   fphi=0.
                   fth=0.
                   qzav=0. & rzav=0. & pzav=0. & uzav=0. & vzav=0.
                   for i=0,nlon-1 do begin
                       st=denth(i,j,l)
                       if st gt 0. then begin
                           ;rwt=0.5*(st+den0) ; density weight
                           rwt=st
                           mpvt=qth(i,j,l)
                           qt=mpvt*thfac
                           ue=uthf(i,j,l)/cosj(j)-u0
                           ve=vthf(i,j,l)/cosj(j)
                           pt=p00*(tth(i,j,l)/thl)^(1./kappa)
                           pe=pt-p0
                           watmp=0.
;
;             Sum up the zonal averages
;
                           qzav=qzav+qt
                           rzav=rzav+st
                           pzav=pzav+pt
                           uzav=uzav+u0+ue
                           vzav=vzav+ve
;
;             Test for Q_k > q at that point
;
                           karr=where(trlevst gt mpvt)
                           k=karr(0)
                           if k gt 0 then begin
                               if emu(k,l) lt 1. then begin
                                   rj=(mpvt-trlevst(k-1)) $
                                     /(trlevst(k)-trlevst(k-1))
                                   mt=(1.-rj)*bsmass(k-1,l)+rj*bsmass(k,l)
                                   ct=(1.-rj)*bscirc(k-1,l)+rj*bscirc(k,l)
                                   cast=(1.-rj)*pcask(k-1,l)+rj*pcask(k,l)
                                   ecast=(1.-rj)*ecask(k-1,l)+rj*ecask(k,l)
                               endif else begin
;               Close to North Pole use B.C. A=M=U=0 at NP.
;
                                   rj=(mpvt-trlevst(kmax)) $
                                     /(qmaxth(l)-trlevst(kmax))
                                   rj=min([rj,1.])
                                   omrj=1.-rj
                                   emut=omrj*emu(kmax,l)+rj*1.
                                   mt=omrj*bsmass(kmax,l)
                                   ucosj=omrj*uthk(kmax,l)*equivcos(kmax,l)
                                   ct=ucosj/(2.*radea)+0.5*omega*(1.-emut*emut)
                                   cast=omrj*pcask(kmax,l) $
                                     +rj*qmaxth(l)*thfac*2.*radea*(-psmass(l))
                                   ecast=ecask(kmax,l) $
                                     +psik(kmax,l)*(mpvt-trlevst(kmax))*thfac
                               endelse
                           endif else begin
                               if k eq 0 then begin
;                                  PV is less than lowest PV level
                                   mt=psmass(l)
                                   ct=pscirc(l)
                                   cast=-2.*radea*pscirc(l)
                                   ecast=ecask(0,l)
                               endif else begin
;                                  PV exceeds values of all PV levels
                                   mt=bsmass(ntrlevst-1,l)
                                   ct=bscirc(ntrlevst-1,l)
                                   cast=pcask(ntrlevst-1,l)
                                   ecast=ecask(ntrlevst-1,l)
                               endelse
                           endelse
                           watmp=c2fix*rwt*dfac*(-qt*(mt-m0)+ct-c0)
;                           waalt=c2fix*st*( -cast $
;                                 -dfac*(c0-qt*(m0-psmass(l))) )
;                           watmp=waalt
                           wa3d(i,j,l)=watmp
                           wazonav=wazonav+watmp
                           walinav=walinav+0.5*den0*den0*dmudq0*(qt-q0)^2
                           wad1=wad1-c2fix*den0*ue*cosj(j) $
                                -c2fix*dfac*(m0-psmass(l))*(st*qt-den0*q0)
                           wagtmp=-(st-den0)*ue*cosj(j)
                           wag=wag+wagtmp
                           h2=thl*cp*((pt/p00)^kappa) $
                             -thl*cp*((p0/p00)^kappa)*(1-kappa+kappa*pt/p0)
                           c2=ecast-ecas0-psi0*(qt-q0)
                           peke=peke+rwt*0.5*(ue*ue+ve*ve)
                           pes=pes+c2fix*rwt*c2
                           peg=peg+(st-den0)*u0*ue
                           ; full APE
;                           pep=pep+rwt*h2+0.5*pe*pe*d2hdpdth/ga
                           ; linearised APE (quadratic form)
                           pep=pep+0.5*pe*pe*cp*((p0/p00)^kappa)*kappa/(ga*p0)
;                           if l eq 55 and j eq 125 then begin
;                               print,i,pep,den0,h2,pe,p0,d2hdpdth
;                           endif
                           ped1=ped1+c2fix*(den0*u0*ue+psi0*(st*qt-den0*q0))
                           fphi=fphi+ve*(watmp+wagtmp) $
                             -den0*ue*ve*cosj(j)
                           fth=fth+pe*dmdlam(i,j,l)/(radea*ga)
;
; Include to treat "interior" on Underworld surface as H_d term.
; Global integral of modified term is equivalent to integral
; of unmodified term plus boundary integral H_b (but less accurate).
;
;                           if jelmax lt nlat then begin
;                               peke=peke+den0*u0*ue $
;                                          +psi0*((st-den0)*qt+den0*(qt-q0))
;                           endif
                        endif
;                  End of longitude loop
                   endfor
;
;                  Now copy zonal integral of wave activity into 2-D arrays.
;                   
;                   if m0*c0 eq 0 then begin
;                       print,'WARNING at location j, l =',j,l
;                       print,m0,c0,ecas0,psi0
;                   endif
                   waden(j,l)=wazonav*dlon
                   wadenlin(j,l)=walinav*dlon
                   wagrav(j,l)=wag*dlon
                   peden(j,l)=peke*dlon
                   pew(j,l)=pes*dlon
                   pegrav(j,l)=peg*dlon
                   peape(j,l)=pep*dlon
                   fluxphi(j,l)=fphi*dlon
                   fluxth(j,l)=fth*dlon
                   if ldterms eq 0 then begin
                       wad(j,l)=wad1*dlon
                       ped(j,l)=ped1*dlon
                   endif
                   if j eq jelmax-1 then begin
                       psi0b=psi0b+psi0
                       m0b=m0b+m0
                   endif
;
;                  Calculate wave activity associated with zonal mean
;
                   if mselect gt 0 then begin
                      qzav=qzav/nlon
                      rzav=rzav/nlon
                      pzav=pzav/nlon
                      uzav=uzav/nlon
                      vzav=vzav/nlon
                      qt=qzav
                      mpvt=qt/thfac
                      st=rzav
                      ;rwt=0.5*(st+den0) ; density weight
                      rwt=st
                      ue=uzav-u0
                      ve=vzav
                      pt=pzav
                      pe=pt-p0
;
;                      for Q_k > [q] at that point
;
                      karr=where(trlevst gt mpvt)
                      k=karr(0)
                      if k gt 0 then begin
                         if emu(k,l) lt 1. then begin
                            rj=(mpvt-trlevst(k-1)) $
                               /(trlevst(k)-trlevst(k-1))
                            mt=(1.-rj)*bsmass(k-1,l)+rj*bsmass(k,l)
                            ct=(1.-rj)*bscirc(k-1,l)+rj*bscirc(k,l)
                            cast=(1.-rj)*pcask(k-1,l)+rj*pcask(k,l)
                            ecast=(1.-rj)*ecask(k-1,l)+rj*ecask(k,l)
                         endif else begin
;               Close to North Pole use B.C. A=M=U=0 at NP.
;
                            rj=(mpvt-trlevst(kmax)) $
                               /(qmaxth(l)-trlevst(kmax))
                            rj=min([rj,1.])
                            omrj=1.-rj
                            emut=omrj*emu(kmax,l)+rj*1.
                            mt=omrj*bsmass(kmax,l)
                            ucosj=omrj*uthk(kmax,l)*equivcos(kmax,l)
                            ct=ucosj/(2.*radea)+0.5*omega*(1.-emut*emut)
                            cast=omrj*pcask(kmax,l) $
                                 +rj*qmaxth(l)*thfac*2.*radea*(-psmass(l))
                            ecast=ecask(kmax,l) $
                                  +psik(kmax,l)*(mpvt-trlevst(kmax))*thfac
                         endelse
                      endif else begin
                         if k eq 0 then begin
;                                  PV is less than lowest PV level
                            mt=psmass(l)
                            ct=pscirc(l)
                            cast=-2.*radea*pscirc(l)
                            ecast=ecask(0,l)
                         endif else begin
;                                  PV exceeds values of all PV levels
                            mt=bsmass(ntrlevst-1,l)
                            ct=bscirc(ntrlevst-1,l)
                            cast=pcask(ntrlevst-1,l)
                            ecast=ecask(ntrlevst-1,l)
                         endelse
                      endelse
                      wazonav=c2fix*rwt*dfac*(-qt*(mt-m0)+ct-c0)
                      walinav=0.5*den0*den0*dmudq0*(qt-q0)^2
                      wad1=-c2fix*den0*ue*cosj(j) $
                           -c2fix*dfac*(m0-psmass(l))*(st*qt-den0*q0)
                      wag=-(st-den0)*ue*cosj(j)
                      h2=thl*cp*((pt/p00)^kappa) $
                         -thl*cp*((p0/p00)^kappa)*(1-kappa+kappa*pt/p0)
                      c2=ecast-ecas0-psi0*(qt-q0)
                      peke=rwt*0.5*(ue*ue+ve*ve)
                      pes=c2fix*rwt*c2
                      peg=(st-den0)*u0*ue
                      pep=rwt*h2+0.5*pe*pe*d2hdpdth/ga
                      ped1=c2fix*(den0*u0*ue+psi0*(st*qt-den0*q0))
                      fphi=ve*(wazonav+wag) $
                           -den0*ue*ve*cosj(j)
                      fth=0. ; cannot be a zonally symmetric part
;
;                  Now take off the wave activity of the zonal mean if
;                  mselect > 0
;
;                      print,'Interior latitude ',j,m,' zonal average WA'
;                      print,wazonav,wag,peke,pes,pep,peg
                      waden(j,l)=waden(j,l)-wazonav*nlon*dlon
                      wadenlin(j,l)=wadenlin(j,l)-walinav*nlon*dlon
                      wagrav(j,l)=wagrav(j,l)-wag*nlon*dlon
                      peden(j,l)=peden(j,l)-peke*nlon*dlon
                      pew(j,l)=pew(j,l)-pes*nlon*dlon
                      pegrav(j,l)=pegrav(j,l)-peg*nlon*dlon
                      peape(j,l)=peape(j,l)-pep*nlon*dlon
                      fluxphi(j,l)=fluxphi(j,l)-fphi*nlon*dlon
                      fluxth(j,l)=fluxth(j,l)-fth*nlon*dlon
                      if ldterms eq 0 then begin
                         wad(j,l)=wad(j,l)-wad1*nlon*dlon
                         ped(j,l)=ped(j,l)-ped1*nlon*dlon
                      endif
                   endif
                   if mselect eq -999 then begin
;
;                     More accurate for Fourier filtered data to 
;                     stick with the linearised expressions for
;                     interior wave activity, since meridional
;                     displacements are small compared with mu_e spacing.
;
                      waden(j,l)=wadenlin(j,l)
                      pew(j,l)=-wadenlin(j,l)*u0/cosj(j)
                   endif

;              End of latitude loop for interior D-bar domain
               endfor
               if emus(l) gt 0. then begin
;
; Find boundary integral contribution.
; jelb=jelmax if waves are large enough for grid points between
; domain boundaries of D-bar and D0. In this case boundary is
; taken to be at midpoint (mu(jelmax)+mu(jelmax-1))/2.
; Otherwise jelb=jelmax-1 and boundary used is mu(jelb).
;
                   if jelb eq jelmax and jelmax gt 0 then begin
                       phib=asin(muh(jelmax))
                       ub0=0.5*(uth0(jelmax-1,l)+uth0(jelmax,l))*cos(phib)
                       ubt=0.5*(uthf(*,jelmax-1,l)+uthf(*,jelmax,l))
                       usum=total(ubt-ub0)*dlon
                       wab(l)=-dfac*(m0b-psmass(l))*usum
                       peb(l)=0.5*psi0b*usum
                   endif else begin
                       if jelb ge 0 then  begin
                           usum=total(uthf(*,jelb,l)-uth0(jelb,l)*cosj(jelb))*dlon
                           wab(l)=-dfac*(m0b-psmass(l))*usum
                           peb(l)=psi0b*usum
                       endif
                   endelse
               endif
           endif
       endfor
;
;      Find bottom boundary pseudoenergy term and
;      boundary term at domain top (approx z0top by z0 on top level).
;
       lowerbound=0
       pe=dblarr(nlon)
       for j=0,nlatwa-1 do begin
           pe(*)=0.
           pszon=p00*total(plb(*,j))/nlon
           psoffset=pszon-ps0(j)
;           psoffset=0.

           case lowerbound of
           0: begin
;
;          Lower boundary defined as an eta-level.
;          For this term choose to define perturbation relative to
;          zonal mean or the MLM background state using psoffset.
;
               pbg=ps0(j)
               pe(*)=plb(*,j)*p00-pbg-psoffset
               dhdpfac=rdgas*ts0(j)*((pbg/p00)^kappa)/(2.*ga*pbg)
           end
           1: begin
;
;          Lower boundary given by D-bar domain. 
;          Calculate pressure perturbation on isentropic surface tsmax.
;
               tsmax=max([thetalb(*,j),ts0(j)])
               thgr=where(thlev gt tsmax)
               lp=thgr(0)
               lm=max([lp-1,0])
               lp=lm+1
               r=(tsmax-thlev(lm))/(thlev(lp)-thlev(lm))
               omr=1.-r
               pm0=omr*pr0(j,lm)+r*pr0(j,lp)
               for i=0,nlon-1 do begin
                   prm=p00*(tth(i,j,*)/thlev(*))^(1./kappa)
                   pm=omr*prm(lm)+r*prm(lp)
                   pe(i)=pm-pm0
               endfor           
               dhdpfac=rdgas*tsmax*((pm0/p00)^kappa)/(2.*ga*pm0)
           end
           2: begin
;
;          Lower boundary given by intersection domain.
;          Note thlev increasing here.
;
               for i=0,nlon-1 do begin
                   if thetalb(i,j) lt ts0(j) then begin
;
;          Lower boundary of intersection domain given by background
;          state. Calculate pressure perturbation on isentropic
;          surface ts0(j).
;
                       thgr=where(thlev gt ts0(j))
                       lp=thgr(0)
                       prm=p00*(tth(i,j,*)/thlev(*))^(1./kappa)
                       lm=max([lp-1,0])
                       lp=lm+1
                       r=(ts0(j)-thlev(lm))/(thlev(lp)-thlev(lm))
                       omr=1.-r
                       pm=omr*prm(lm)+r*prm(lp)
                       pe(i)=pm-ps0(j)-psoffset
                   endif else begin
;
;          Lower boundary of intersection domain given by 3D state.
;          Calculate pressure perturbation on isentropic surface thetalb(i,j).
;
                       thgr=where(thlev gt thetalb(i,j))
                       lp=thgr(0)
                       lm=max([lp-1,0])
                       lp=lm+1
                       r=(thetalb(i,j)-thlev(lm))/(thlev(lp)-thlev(lm))
                       omr=1.-r
                       pm0=omr*pr0(j,lm)+r*pr0(j,lp)
                       pe(i)=plb(i,j)*p00-pm0-psoffset
                   endelse
               endfor
               dhdpfac=rdgas*ts0(j)*((ps0(j)/p00)^kappa)/(2.*ga*ps0(j))
           end
           3: begin
;
;          Alternatively, define pressure perturbation by referencing
;          to background state position of theta_s contour, allowing
;          for its meridional displacement along lower boundary.
;          Note that j is ordered NP to equator.
;          If thetalb is out of range of ts0, then use ends. 
;
              for i=0,nlon-1 do begin
                 tsfull=thetalb(i,j)
                 latgr=where(ts0 gt tsfull)
                 jp=latgr(0)
                 if jp eq -1 then jp=nlat-1
                 jm=max([jp-1,0])
                 jp=jm+1
                 r=(tsfull-ts0(jm))/(ts0(jp)-ts0(jm))
                 omr=1.-r
                 pms0=omr*ps0(jm)+r*ps0(jp)
                 pe(i)=plb(i,j)*p00-pms0
              endfor
              dhdpfac=rdgas*ts0(j)*((ps0(j)/p00)^kappa)/(2.*ga*ps0(j))
           end
           endcase

           pezon=total(pe)/nlon
           pe=pe-pezon
;           bot=-pe*zs0(j)/ga+pe*pe*dhdpfac
           bot=pe*pe*dhdpfac
           pet(j)=total(bot)*dlon

           ptbav=pzon(j)
           pe=p00*ptop(*,j)-ptbav
           temp0top=thtop*((ptbav/p00)^kappa)
;           zt0=(gom0(j,nthlev-1)-cp*temp0(j,nthlev-1))/ga
           dhdpfac=rdgas*temp0top/(2.*ga*ptbav)
;           top=pe*zt0-pe*pe*dhdpfac
           top=-pe*pe*dhdpfac
           petop(j)=total(top)*dlon
       endfor
;
;      End of wave activity density calculations (and loop over theta_m).
;
;      Calculate pseudomomentum flux divergence.
;
       for m=1,nthlim-2 do begin
          for j=1,nlatwa-2 do begin
             fluxdiv(j,m)=(1./(radea*cosj(j))) $
                         *(fluxphi(j+1,m)*cosj(j+1)-fluxphi(j-1,m)*cosj(j-1))$
                         /((latwa(j+1)-latwa(j-1))*!dpi/180.)$
                         +(fluxth(j,m+1)-fluxth(j,m-1))/(thlev(m+1)-thlev(m-1))
          endfor
       endfor
;
;      Find contributions to global integral wave activity.
;
       masstot=0.
       stratmass=0.
       checkmass1=0.
       checkmass2=0.
       checkmass3=0.
       bsmnew=dblarr(ntrlevst,nthlev)
       bsmnew(*,*)=0.
       psmnew=dblarr(nthlev)
       psmnew(*)=0.
       pscheck=dblarr(nlatwa)
       pscheck(*)=0.
       watot=0.
       wagtot=0.
       wadtot=0.
       waetot=0.
       wabtot=0.
       wae2tot=0.
       pee2tot=0.
       petot=0.
       pewtot=0.
       pegtot=0.
       peapetot=0.
       pedtot=0.
       peetot=0.
       pettot=0.
       petoptot=0.
       pebtot=0.
       ; nonz=where(emus lt xtropic)
       ; print,'emus away from equator = ',emus(nonz(0)-1),emus(nonz(0:1))
       ; xnext=emus(nonz(0)-1)
       xnext=xtropic
       for l=0,nthlev-1 do begin
           thl=thlev(l)
           dth=thlevh(l+1)-thlevh(l)
           for j=0,nlatwa-1 do begin
               weight=radea*radea*dmu(j)*dth
               watot=watot+waden(j,l)*weight
               wagtot=wagtot+wagrav(j,l)*weight
               wadtot=wadtot+wad(j,l)*weight
               waetot=waetot+wae(j,l)*weight
               wae2tot=wae2tot+waetmp(j,l)*weight
               pee2tot=pee2tot+peetmp(j,l)*weight
               petot=petot+peden(j,l)*weight
               pewtot=pewtot+pew(j,l)*weight
               pegtot=pegtot+pegrav(j,l)*weight
               peapetot=peapetot+peape(j,l)*weight
               pedtot=pedtot+ped(j,l)*weight
               peetot=peetot+pee(j,l)*weight
               checkmass1=checkmass1+total(denth(*,j,l))*weight*dlon
               checkmass2=checkmass2+sigma0(j,l)*weight*2.*!dpi
               pscheck(j)=pscheck(j)+sigma0(j,l)*dth*ga
           endfor
           wabtot=wabtot+wab(l)*radea*dth
           pebtot=pebtot+peb(l)*radea*dth
           interior=where(ts0 lt thl) 
           if interior(0) gt -1 then begin
               lmass=sigma0(*,l)*dmu/2.         
               psmnew(l)=total(lmass(interior))
               masstot=masstot+psmass(l)*4.*!dpi*radea*radea*dth
               checkmass3=checkmass3+psmnew(l)*4.*!dpi*radea*radea*dth
           endif

           qmax=max(pv0(*,l),min=qmin)
           qsmall=where(trlevst le qmin)
           nsmall=n_elements(qsmall)
           if qsmall(0) eq -1 then begin
               qsmall=0
           endif
           bsmnew(qsmall,l)=psmnew(l)
           for k=nsmall,ntrlevst-1 do begin
               mpv=trlevst(k)
               highq=where(pv0(*,l) gt mpv)
               if highq(0) gt -1 then begin
                   lmass=sigma0(*,l)*dmu/2.
                   bsmnew(k,l)=total(lmass(highq))
               endif           
           endfor
       endfor
       pettot=radea*radea*total(pet*dmu)
       petoptot=radea*radea*total(petop*dmu)
       stratmass=total(pzon*dmu(0:nlatwa-1))
       pstot=total(pszav(0:nlatwa-1)*dmu(0:nlatwa-1))
       if lbound eq 1 then begin
           pt=alev(nlev-nbound)*pstot*p00
       endif else begin
           pt=p00
       endelse
       massnh=(pt-stratmass)*2*!dpi*radea*radea/ga

       print,' '
       print,'Mass check as fraction of correct NH mass below domain top'
       print,'from psmass, denth, sigma0, psmnew arrays'
       print,masstot/massnh,checkmass1/massnh $
            ,checkmass2/massnh,checkmass3/massnh

       print,' '
       print,'Global integral pseudomomentum terms'
       print,'WA interior / hemispheric mass = ',watot/masstot
       print,'grav interior / WAint  = ',wagtot/watot
       print,'d-term / WAint        = ',wadtot/watot
       print,'e-term / WAint        = ',waetot/watot
       print,'boundary-term / WAint  = ',wabtot/watot
       print,'Sum / WAint            = ' $
            ,(watot+wagtot+wadtot+waetot+wabtot)/watot
       print,' '
       print,'Global integral pseudoenergy terms'
       print,'KE interior / hemispheric mass = ',petot/masstot
       print,'wind-weighted interior / KEint  = ',pewtot/petot
       print,'grav interior / KEint  = ',pegtot/petot
       print,'APE / KEint            = ',peapetot/petot
       print,'d-term / KEint        = ',pedtot/petot
       print,'e-term / KEint        = ',peetot/petot
       print,'e1-term / KEint        = ',(peetot-pee2tot)/petot
       print,'e2-term / KEint        = ',pee2tot/petot
       print,'boundary-term / KEint  = ',pebtot/petot
       print,'ground-term / KEint    = ',pettot/petot
       print,'top-term / KEint       = ',petoptot/petot
       print,'Sum / KEint            = ' $
     ,(petot+pewtot+pegtot+peapetot+pedtot+peetot+pebtot+pettot+petoptot)/petot
    endif
;
;  End of wave activity calculations
;
;  Plot PV calculated by two different methods.
;
   if lpvplot ge 1 then begin
       if nlat eq nlatnh then begin
           projchoice=3
       endif else begin
           projchoice=1
       endelse
       etachoice=thchoice
       print,'PV from eta-calc, PV from isentropic vorticity and difference'
       print,'Maximum values (PVU) = ',pvu*max(pvth(*,*,etachoice)) $
         ,pvu*max(qth(*,*,etachoice)) $
         ,pvu*max(pvth(*,*,etachoice)-qth(*,*,etachoice))
       print,'Minimum values (PVU) = ',pvu*min(pvth(*,*,etachoice)) $
         ,pvu*min(qth(*,*,etachoice)) $
         ,pvu*min(pvth(*,*,etachoice)-qth(*,*,etachoice))
       attrchoice=ipv
       lsurf=0 & lisen=1
       
       pvname='xyp'+strtrim(projchoice,2)+'_l'+strtrim(etachoice,2) $
         +'_a'+strtrim(attrchoice,2)+'_'+strtrim(idatadate,2)
       erase
       
       pvvar=dblarr(nlon,nlatwa,nthlev)
       if lmode gt 1 then begin
           for m=0,nthlev-1 do begin
               for j=0,nlatwa-1 do begin
;                   pvvar(*,j,m)=bsmodpv(j,m)
                   pvvar(*,j,m)=qth(*,j,m)
;                   pvvar(*,j,m)=0.01*sigma0(j,m)*qth(*,j,m)
;                   pvvar(*,j,m)=0.01*sigma0(j,m)*(qth(*,j,m)-pv0(j,m))
;                   pvvar(*,j,m)=0.01*(denth(*,j,m)*qth(*,j,m) $
;                                     -sigma0(j,m)*pv0(j,m))
               endfor
           endfor
           pvvar=pvu*pvvar
           units='PV (PVU)'
;           units='PV anomaly (PVU)'
;           units='PV*r0  (10!e-4!ns!e-1!n)'
;           units='PV anomaly (10!e-4!ns!e-1!n)'
;           units='Vorticity anomaly (10!e-4!ns!e-1!n)'

           case lpvplot of
               0: begin

               end
               1: begin
                   lev0=-2.00 & dlev=0.5 & llog=0 & lfill=0 & nozer=0
;                   lev0=-1.00 & dlev=0.10 & llog=0 & lfill=2 & nozer=0
                   plotxy,idatadate,pvvar,longitude,latwa,alev $
                     ,lsurf,lisen,thlev,attrchoice,units $
                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
;                   erase
;                   attrchoice=itheta
;                   units='T  (K)'
;                   lev0=200. & dlev=5. & llog=0 & lfill=2 & nozer=0
;                   plotxy,idatadate,tth(*,0:nlatwa-1,*),longitude,latwa,alev $
;                     ,lsurf,lisen,thlev,attrchoice,units $
;                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
               end
               2: begin
                   attrchoice=ipv+1
                   units='PV diff (PVU)'
                   lev0=-0.75 & dlev=0.05
                   plotxy,idatadate,pvu*(qth-pvth),longitude,latitude,alev $
                     ,lsurf,lisen,thlev,attrchoice,units $
                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
               end
               4: begin
                   attrchoice=ipv+1
                   units='wave activity density'
                   pvvar=wa3d
                   lev0=min(pvvar(*,*,etachoice),max=maxwa) & dlev=(maxwa-lev0)/nclev
                   print,'lev0 dlev = ',lev0,dlev
                   plotxy,idatadate,wa3d,longitude,latitude,alev $
                     ,lsurf,lisen,thlev,attrchoice,units $
                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
;                   plotxy,idatadate,wa3dqm0,longitude,latitude,alev $
;                     ,lsurf,lisen,thlev,attrchoice,units $
;                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
;                   pvvar=wa3dct
;                   lev0=min(pvvar(*,*,etachoice),max=maxwa) & dlev=(maxwa-lev0)/nclev
;                   plotxy,idatadate,wa3dct,longitude,latitude,alev $
;                     ,lsurf,lisen,thlev,attrchoice,units $
;                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
;                   plotxy,idatadate,wa3dc0,longitude,latitude,alev $
;                     ,lsurf,lisen,thlev,attrchoice,units $
;                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
               end
               5: begin
                   attrchoice=iu
                   units='u cos(lat)  (m/s)'
                   lev0=-32. & dlev=2. & llog=0 & lfill=2 & nozer=0
                   plotxy,idatadate,uthf(*,0:nlatwa-1,*),longitude,latwa,alev $
                     ,lsurf,lisen,thlev,attrchoice,units $
                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
               end
               6: begin
                   attrchoice=iv
                   units='v cos(lat)  (m/s)'
                   lev0=-16. & dlev=1. & llog=0 & lfill=2 & nozer=0
                   plotxy,idatadate,vthf(*,0:nlatwa-1,*),longitude,latwa,alev $
                     ,lsurf,lisen,thlev,attrchoice,units $
                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
               end
               12: begin
                   attrchoice=iv
                   units='streamfunction  (m2/s)'
                   lev0=minmont & dlev=dmont & llog=0 & lfill=2 & nozer=0
;                   plotxy,idatadate,stream,longitude,latwa,alev $
;                     ,lsurf,lisen,thlev,attrchoice,units $
;                     ,lev0,dlev,nclev,llog,lfill,nozer,etachoice,projchoice
               end
           endcase
           elats=asin(emus(etachoice))*180./!pi
           oplot,[0,359],[elats,elats],thick=4
           
           if hard eq 1 then begin
               device,/close
               print,' '
               print,' WARNING: enter .cont to output waden and other plots'
               print,'*****************************************************'
               stop
           endif
           
           if hard eq 2 then begin
               tvlct, red, green, blue, /get
               case visname of
                   'PseudoColor': begin
                       scrimage=tvrd()
                   end
                   'TrueColor': begin
                       scrimage=tvrd(0,0,true=1)
                   end
               endcase
;          scrimage=reverse(tvrd(),2)
               pngname=pvname+'.png'
               write_png,pngname,scrimage,red,green,blue
           endif
       endif
   endif
;
; Calculate meridional wind in latitude band on level thchoice
;
   nband=jbands-jbandn+1
   for m=0,nthlim-1 do begin
       j=jbandn
       vband(*,m)=vthf(*,j,m)/cosj(j)
       for j=jbandn+1,jbands do begin
           vband(*,m)=vband(*,m)+vthf(*,j,m)/cosj(j)
       endfor
   endfor
   vband=vband/nband
;
; Plots
;
   if hard eq 1 and lmode gt 1 then begin
      psname='bsout_'+fstembi+strtrim(idatadate,2)
      erase
      if hard eq 1 then begin
          set_plot,'ps'
          device, /landscape, /color, bits_per_pixel=8
          device, file=psname+'.ps'
;          loadct,ctab
      endif
      !p.multi=[0,3,2]
      !y.omargin=[0,4]
      !p.charsize=1.5

;      topth=max(thlev)
      topth=400.
      trlev=trlevst*pvu
;      contour,bsmass(*,*),trlev,thlev $
;             ,title='Mass enclosed by PV contours ' $
;             ,xtitle='Modified PV  (PVU)' $
;             ,xrange=[trlev(1),max(trlev)],xstyle=1,/xlog $
;             ,ytitle='Potential temperature (K)' $
;             ,yrange=[thlev0,max(thlev)],ystyle=1 $
;             ,levels=clevs,c_linestyle=clines

      cmax=max(bscirc,min=cmin)
      cilev=cmin+((cmax-cmin)/nclev)*findgen(nclev)
;      contour,bscirc(*,*),trlev,thlev $
;             ,title='Circulation integral ' $
;             ,xtitle='Modified PV  (PVU)' $
;             ,xrange=[trlev(1),max(trlev)],xstyle=1,/xlog $
;             ,ytitle='Potential temperature (K)' $
;             ,yrange=[thlev0,max(thlev)],ystyle=1 $
;             ,levels=cilev

      clabs=trlev
      clabs(*)=1
      contour,bsmodpv(*,*)*pvu,latwa,thlev $
             ,title='Background state modified PV  (PVU)' $
             ,xtitle='Latitude' $
             ,xrange=[0.,90.],xstyle=1 $
             ,ytitle='Potential temperature (K)' $
;             ,ytitle='Height  (km)' $
             ,yrange=[thlev0,topth],ystyle=1 $
             ,levels=trlev,c_labels=clabs

      ulevs=-40.+10.*findgen(16)
      ulabs=ulevs
      ulabs(*)=1
      uline=(ulevs lt 0)
      contour,uthinv(*,*),rlatinv,thlev $
             ,title='Inverter-grid u  (m s!e-1!n)' $
             ,xtitle='Latitude' $
             ,xrange=[0.,90.],xstyle=1 $
             ,ytitle='Potential temperature (K)' $
;             ,ytitle='Height  (km)' $
             ,yrange=[thlev0,topth],ystyle=1 $
             ,levels=ulevs,c_labels=ulabs,c_linestyle=uline

      slevs=20.*findgen(16)
      slabs=slevs
      slabs(*)=1
      sline=(slevs lt 0)
      contour,sigmainv(*,*),rlatinv,thlev $
             ,title='Inverter-grid density  (kg m!e-2!n K!e-1!n)' $
             ,xtitle='Latitude' $
             ,xrange=[0.,90.],xstyle=1 $
             ,ytitle='Potential temperature (K)' $
;             ,ytitle='Height  (km)' $
             ,yrange=[thlev0,topth],ystyle=1 $
             ,levels=slevs,c_labels=slabs,c_linestyle=sline

;      contour,sigma0(*,*),latitude,thlev $
;             ,title='Background state density  (kg m!e-2!n K!e-1!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,max(thlev)],ystyle=1 $
;             ,levels=slevs,c_labels=slabs,c_linestyle=sline

      cmax=max(gominv,min=cmin)
      glevs=cmin+((cmax-cmin)/nclev)*findgen(nclev)
;      glevs=20.*findgen(16)
      glabs=glevs
      glabs(*)=1
      gline=(glevs lt 0)
;      contour,gominv(*,*),rlatinv,thlev $
;             ,title='Inverter-grid Montgomery  (m!e2!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,max(thlev)],ystyle=1 $
;             ,levels=glevs,c_labels=glabs,c_linestyle=gline

      cmax=max(prinv,min=cmin)
      plevs=cmin+((cmax-cmin)/nclev)*findgen(nclev)
;      plevs=20.*findgen(16)
      plabs=plevs
      plabs(*)=1
      pline=(plevs lt 0)
 ;     contour,prinv(*,*),rlatinv,thlev $
 ;            ,title='Inverter-grid pressure  (kg m!e-1!n s!e-2!n)' $
 ;            ,xtitle='Latitude' $
 ;            ,xrange=[0.,90.],xstyle=1 $
 ;            ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
 ;            ,yrange=[thlev0,max(thlev)],ystyle=1 $
 ;            ,levels=plevs,c_labels=plabs,c_linestyle=pline

      cmax=max(psik,min=cmin)
      cilev=cmin+((cmax-cmin)/nclev)*findgen(nclev)
;      contour,psik(*,*),trlev,thlev $
;             ,title='Streamfunction ' $
;             ,xtitle='Modified PV  (PVU)' $
;             ,xrange=[trlev(1),max(trlev)],xstyle=1,/xlog $
;             ,ytitle='Potential temperature (K)' $
;             ,yrange=[thlev0,max(thlev)],ystyle=1 $
;             ,levels=cilev

      cmax=max(ecask,min=cmin)
      cilev=cmin+((cmax-cmin)/nclev)*findgen(nclev)
;      contour,ecask(*,*),trlev,thlev $
;             ,title='Casimir ' $
;             ,xtitle='Modified PV  (PVU)' $
;             ,xrange=[trlev(1),max(trlev)],xstyle=1,/xlog $
;             ,ytitle='Potential temperature (K)' $
;             ,yrange=[thlev0,max(thlev)],ystyle=1 $
;             ,levels=cilev

      wamax=max(waden)
      walevs=(-10+findgen(41))*400
      walabs=walevs
      walabs(*)=1
      waline=(walevs lt 0)
      contour,waden(*,*),latwa,thlev $
             ,title='nonlinear WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
             ,xtitle='Latitude' $
             ,xrange=[0.,90.],xstyle=1 $
             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
             ,yrange=[thlev0,topth],ystyle=1 $
             ,levels=walevs,c_labels=walabs,c_linestyle=waline

;      contour,wagrav(*,*),latwa,thlev $
;             ,title='g-term WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=walevs,c_labels=walabs,c_linestyle=waline

;      contour,wad(*,*),latwa,thlev $
;             ,title='d-term WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=walevs,c_labels=walabs,c_linestyle=waline

;      contour,waetmp(*,*),latwa,thlev $
;             ,title='d2-term WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=walevs,c_labels=walabs,c_linestyle=waline

;      contour,wae(*,*),latwa,thlev $
;             ,title='e-term WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=walevs,c_labels=walabs,c_linestyle=waline

;      contour,peetmp(*,*),latwa,thlev $
;             ,title='e2-term WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=walevs,c_labels=walabs,c_linestyle=waline

      contour,wadenlin(*,*),latwa,thlev $
             ,title='linear WA density  (kg m!e-1!n K!e-1!n s!e-1!n)' $
             ,xtitle='Latitude' $
             ,xrange=[0.,90.],xstyle=1 $
             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
             ,yrange=[thlev0,topth],ystyle=1 $
             ,levels=walevs,c_labels=walabs,c_linestyle=waline

;      l=20
;      plot,latwa,waden(*,l) $
;        ,title='WA density on '+strtrim(floor(thlev(l)),2)+'K' $
;        ,xtitle='Latitude',xrange=[0.,90.],xstyle=1 $
;        ,ytitle='WA density (kg m!e-1!n K!e-1!n s!e-1!n)',yrange=[min(walevs),max(walevs)]
;      oplot,latitude,wagrav(*,l)
;      oplot,latitude,wad(*,l),linestyle=1
;      oplot,latitude,waetmp(*,l),linestyle=1
;      oplot,latitude,wae(*,l)
;      oplot,latitude,peetmp(*,l)

      pemax=max([peden,-pew])
      pelevs=pemax*(-1.+(0.5+findgen(21))/10.)
      pelabs=pelevs
      pelabs(*)=1
      peline=(pelevs lt 0)
;      contour,peden(*,*),latwa,thlev $
;             ,title='KE density  (kg K!e-1!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=pelevs,c_labels=pelabs,c_linestyle=peline

;      contour,-pew(*,*),latwa,thlev $
;             ,title='wind-term PE density  (kg K!e-1!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=pelevs,c_labels=pelabs,c_linestyle=peline

;      contour,pegrav(*,*),latwa,thlev $
;             ,title='g-term PE density  (kg K!e-1!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=pelevs,c_labels=pelabs,c_linestyle=peline

;      contour,peape(*,*),latwa,thlev $
;             ,title='APE density  (kg K!e-1!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=pelevs,c_labels=pelabs,c_linestyle=peline

;      contour,ped(*,*),latwa,thlev $
;             ,title='d-term PE density  (kg K!e-1!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=pelevs,c_labels=pelabs,c_linestyle=peline

;      contour,pee(*,*),latwa,thlev $
;             ,title='e-term PE density  (kg K!e-1!n s!e-2!n)' $
;             ,xtitle='Latitude' $
;             ,xrange=[0.,90.],xstyle=1 $
;             ,ytitle='Potential temperature (K)' $
;;             ,ytitle='Height  (km)' $
;             ,yrange=[thlev0,topth],ystyle=1 $
;             ,levels=pelevs,c_labels=pelabs,c_linestyle=peline

      xyouts,0.4,0.96,strtrim(idatadate,2),/normal,charsize=1.5

      if hard eq 2 then begin
          scrimage=reverse(tvrd(),2)
          write_tiff,psname+'.tiff',scrimage,red=red,green=green,blue=blue
      endif
      if hard eq 1 then begin
          device,/close
      endif
   endif

   if lmode gt 1 then begin
;       print,'Stopped just before reversing arrays for output'
;       print,'Re-ordering theta so that it runs from top to bottom'
;       print,'WARNING: have not output results!' & stop
       lboundout=1 ; select 1 to output lower boundary of bg state
                   ; Note: this was only added on 30/7/2015
       pv0=reverse(pv0,2)
       qbaryj=reverse(qbaryj,2)
       uth0=reverse(uth0,2)
;       uth0=reverse(bsucosj,2)
       sigma0=reverse(sigma0,2)
       pr0=reverse(pr0,2)
       zavucosj=reverse(zavucosj,2)
       waden=reverse(waden,2) & wagrav=reverse(wagrav,2)
       wad=reverse(wad,2) & wae=reverse(wae,2) & wab=reverse(wab)
       peden=reverse(peden,2) & pegrav=reverse(pegrav,2) & pew=reverse(pew,2)
       ped=reverse(ped,2) & pee=reverse(pee,2) & peb=reverse(peb)
       peape=reverse(peape,2) 
;       pet=pet+petop
       pet=reverse(pet)
       fluxphi=reverse(fluxphi,2) & fluxth=reverse(fluxth,2) 
       fluxdiv=reverse(fluxdiv,2)
       vband=reverse(vband,2)

       if mselect eq 0 then begin
           wafile=fpatho+'waden_'+fstembi+smode+strtrim(idatadate,2)
       endif else begin
           wafile=fpatho+'waden_'+strtrim(expstem,2) $
                 +'_m'+strtrim(mselect,2)+'_'+strtrim(idatadate,2)
       endelse
       openw,2,wafile
       printf,2,' DATA BASE TIME IS ',idatadate,format='(A19,I10)'
       printf,2,nlatwa,' LATITUDES ON GAUSSIAN GRID AT',format='(I4,A)'
       printf,2,latwa,format='(8E13.5)'
       printf,2,nthlev,' ISENTROPIC LEVELS AT',format='(I4,A)'
;       printf,2,thlevrev,format='(10F8.3)'
;       printf,2,thtop,' IS TOP BOUNDARY IN ISENTROPIC COORDS',format='(F8.3,A)'
       printf,2,thlevrev,format='(10F8.2)'
       printf,2,thtop,' IS TOP BOUNDARY IN ISENTROPIC COORDS',format='(F8.2,A)'
       printf,2,thsmax,thsmin,' MAX AND MIN VALUES OF SURFACE THETA' $
         ,format='(2E13.5,A)'
       if lboundout eq 1 then begin
          printf,2,' Z ALONG LOWER BOUNDARY OF BACKGROUND'   
          printf,2,zs0,format='(8E13.5)'
          printf,2,' PRESSURE ALONG LOWER BOUNDARY OF BACKGROUND'   
          printf,2,ps0,format='(8E13.5)'
          printf,2,' THETA ALONG LOWER BOUNDARY OF BACKGROUND'   
          printf,2,ts0,format='(8E13.5)'
          printf,2,' ZONAL WIND ALONG LOWER BOUNDARY OF BACKGROUND'   
          printf,2,us0,format='(8E13.5)'
       endif
       printf,2,' BACKGROUND PV IN LATITUDE-THETA COORDINATES'   
       printf,2,pv0,format='(8E13.5)'
       printf,2,' BACKGROUND DENSITY IN LATITUDE-THETA COORDINATES'   
       printf,2,sigma0,format='(8E13.5)'
       printf,2,' BACKGROUND PRESSURE IN LATITUDE-THETA COORDINATES'   
       printf,2,pr0,format='(8E13.5)'
       printf,2,' BACKGROUND PV GRADIENT IN LATITUDE-THETA COORDINATES'   
       printf,2,qbaryj,format='(8E13.5)'
       printf,2,' BACKGROUND u IN LATITUDE-THETA COORDINATES' 
       printf,2,uth0,format='(8E13.5)'
       printf,2,' ZONAL AVE u cos(lat) IN 1/2-LATITUDE-THETA COORDINATES' 
       printf,2,zavucosj,format='(8E13.5)'
       printf,2,' GLOBAL INTEGRALS OF MASS FROM 3D AND INVERTER STATES'
       printf,2,masstot,checkmass2,format='(8E13.5)'
       printf,2,' GLOBALLY AVERAGED PSEUDOMOMENTUM TERMS'
       printf,2,watot,wagtot,wadtot,waetot,wabtot,wae2tot,format='(8E13.5)'
       printf,2,' GLOBALLY AVERAGED PSEUDOENERGY TERMS'
       printf,2,petot,pewtot,pegtot,pedtot,peetot,pebtot,peapetot $
               ,pettot,petoptot,format='(9E13.5)'
       printf,2,$
       ' ZONALLY AVERAGED PSEUDOMOMENTUM DENSITY IN LATITUDE-THETA COORDINATES' 
       printf,2,waden,format='(8E13.5)'
       printf,2,wagrav,format='(8E13.5)'
       printf,2,wad,format='(8E13.5)'
       printf,2,wae,format='(8E13.5)'
       printf,2,wab,format='(8E13.5)'
       printf,2,$
       ' ZONALLY AVERAGED PSEUDOENERGY DENSITY IN LATITUDE-THETA COORDINATES' 
       printf,2,peden,format='(8E13.5)'
       printf,2,pew,format='(8E13.5)'
       printf,2,pegrav,format='(8E13.5)'
       printf,2,ped,format='(8E13.5)'
       printf,2,pee,format='(8E13.5)'
       printf,2,peb,format='(8E13.5)'
       printf,2,peape,format='(8E13.5)'
       printf,2,pet,format='(8E13.5)'
       printf,2,$
       ' ZONALLY AVERAGED PSEUDOMOMENTUM FLUX IN LATITUDE-THETA COORDINATES' 
       printf,2,fluxphi,format='(8E13.5)'
       printf,2,fluxth,format='(8E13.5)'
       printf,2,fluxdiv,format='(8E13.5)'
       printf,2,' MERIDIONAL WIND ALONG LATITUDE BAND'
       printf,2,phibands,phibandn,format='(8E13.5)'
       printf,2,vband,format='(8E13.5)'
       close,2

       if mselect gt 0 then begin
;           print,'Re-ordering theta so that it runs from top to bottom'
           qifour=reverse(qifour,2)
           pifour=reverse(pifour,2)
           uifour=reverse(uifour,2)
           vifour=reverse(vifour,2)
           
           openw,12,fpatho+'fourier_'+strtrim(expstem,2) $
                   +'_m'+strtrim(mselect,2)+'_'+strtrim(idatadate,2)
           printf,12,' DATA BASE TIME IS ',idatadate,format='(A19,I10)'
           printf,12,nlatwa,' LATITUDES ON GAUSSIAN GRID AT',format='(I4,A)'
           printf,12,latwa,format='(8E13.5)'
           printf,12,nthlev,' ISENTROPIC LEVELS AT',format='(I4,A)'
           printf,12,thlevrev,format='(10F8.2)'
           printf,12,thtop,' IS TOP BOUNDARY IN ISENTROPIC COORDS',format='(F8.2,A)'
           printf,12,thsmax,thsmin,' MAX AND MIN VALUES OF SURFACE THETA' $
             ,format='(2E13.5,A)'
           printf,12,' FOURIER COEFFICIENTS FOR ZONAL WAVENUMBER ' $
             ,mselect*moct,format='(A,I4)'           
           printf,12,' COEFFICIENTS FOR SURFACE z, p, theta, u, v'
           printf,12,' STORED ON ANALYSIS LATITUDE POINTS '
           printf,12,zsfour,format='(8E13.5)'
           printf,12,psfour,format='(8E13.5)'
           printf,12,tsfour,format='(8E13.5)'
           printf,12,usfour,format='(8E13.5)'
           printf,12,vsfour,format='(8E13.5)'
           printf,12,' COEFFICIENTS FOR INTERIOR q, p, u, v'
           printf,12,' STORED ON ANALYSIS LATITUDE x THETA POINTS '
           printf,12,qifour,format='(8E13.5)'
           printf,12,pifour,format='(8E13.5)'
           printf,12,uifour,format='(8E13.5)'
           printf,12,vifour,format='(8E13.5)'
           close,12
       endif
   endif

skip:
   if hard eq 1 then begin
       device,/close
   endif   
   thetam=theta
   theta=reform(datai[*,*,*,itheta],nlon,nlat,nlev)
   pv=reform(datai[*,*,*,ipv],nlon,nlat,nlev)
   pv=pv/pvu
   u=reform(datai[*,*,*,iu],nlon,nlat,nlev)
   v=reform(datai[*,*,*,iv],nlon,nlat,nlev)
   psurf=psurfp
   if iheat gt -1 then begin
       thetaadv=reform(datai[*,*,*,iheat],nlon,nlat,nlev)
   endif
   idatadate=iplusdate
   iday=(idatadate-10000*(idatadate/10000))/100
   print,iday


print,'REACHED END SUCCESSFULLY'

end
