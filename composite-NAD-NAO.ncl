
ad "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  f1=addfile("/mnt/d/工作文档/田琦数据/slp.mon.mean3.nc", "r")
    
  wks = gsn_open_wks("png","/mnt/d/工作文档/田琦数据/Comp-NAD-NAO-Patterns") 
  gsn_define_colormap(wks, "cmp_b2r")

  year = ispan(1979,2013,1)
  beginyrmon2=197901;     
  lastyrmon2=201412;         
  time2=f1->time   
  YYYYMM2=cd_calendar(time2,-1)    
  
  rec_beginyrmon2=ind(beginyrmon2.eq.YYYYMM2)    
  rec_lastyrmon2=ind(lastyrmon2.eq.YYYYMM2)      

  slpt=f1->slp(rec_beginyrmon2:rec_lastyrmon2,:,:) 
  t=slpt(2,2,143)
  print(t)
  lat1=(f1->lat)
  lat=-lat1
  lon=(f1->lon)
  timet=rec_beginyrmon2-rec_lastyrmon2
  latinter=fspan(-90.,90.,73)
  loninter=fspan(0.,357.5,144)
  loninter@units="degrees_east"
latinter@units="degrees_north"

  slp1=linint2_Wrap(lon,lat,slpt,True,loninter,latinter,0)

 ; poisson_grid_fill(slp1,True,1,1500,1e-2,0.6,0)



 
  yyyymm2=cd_calendar(slp1&time, -1)
  slp1_lat = f1->lat 

  begincl2=198101
  lastcl2=201012
  iClmStrt2=ind(yyyymm2.eq.begincl2)
  iClmLast2=ind(yyyymm2.eq.lastcl2)

  slpClm = clmMonTLL(slp1(iClmStrt2:iClmLast2,:,:))
  slpAnom=calcMonAnomTLL(slp1,slpClm)

  slp_DJF      = new((/35,73,144/),float)
  slp_DJF!0    = "time"
  slp_DJF&time = year 
  slp_DJF!1    = "lat"
  slp_DJF&lat  = latinter
  slp_DJF!2    = "lon"
  slp_DJF&lon  = loninter

  do i = 0,34
        slp_DJF(i,:,:) = (slpAnom(11+12*i,:,:)+slpAnom(12+12*i,:,:)+slpAnom(13+12*i,:,:))/3
  end do

  f2=addfile("/mnt/d/工作文档/田琦数据/NAD_index.nc", "r")
  NAD=f2->NAD
  NAD_DJF=new(35, float)
  NAD_DJF!0    = "time"
  NAD_DJF&time = year 
  do i = 0,34
        NAD_DJF(i) = (NAD(11+12*i)+NAD(12+12*i)+NAD(13+12*i))/3
  end do

  NAO=asciiread("/mnt/d/工作文档/田琦数据/NAO-CPC79-14.txt", -1, "float")
  NAO_DJF=new(35, float)
  NAO_DJF!0    = "time"
  NAO_DJF&time = year 
  do i = 0,34
        NAO_DJF(i) = (NAO(11+12*i)+NAO(12+12*i)+NAO(13+12*i))/3
  end do

  ;***************************Composite******************************
  ;*******************************************************************
  NAD_NAO_DJF_gt5 = ind(NAD_DJF.gt.0.5.and.NAO_DJF.gt.0.5)
  nnumb=dimsizes(NAD_NAO_DJF_gt5)   

  NAD_NAO_DJF_lt5 = ind(NAD_DJF.lt.-0.5.and.NAO_DJF.lt.-0.5)
  nnumb2=dimsizes(NAD_NAO_DJF_lt5) 

  NAO_DJF_gt5_NAD_lt0 = ind(NAD_DJF.lt.0.and.NAO_DJF.gt.0.5)
  nnumb3=dimsizes(NAO_DJF_gt5_NAD_lt0) 

  NAO_DJF_gt5_NAD_lt5 = ind(NAD_DJF.lt.0.5.and.NAO_DJF.gt.0.5)
  nnumb4=dimsizes(NAO_DJF_gt5_NAD_lt5) 

  NAO_DJF_gt0_NAD_lt0 = ind(NAD_DJF.lt.0.and.NAO_DJF.gt.0)
  nnumb5=dimsizes(NAO_DJF_gt0_NAD_lt0) 
  print(NAO_DJF_gt0_NAD_lt0)

  NAO_DJF_lt0_NAD_gt0 = ind(NAD_DJF.gt.0.and.NAO_DJF.lt.0)
  nnumb6=dimsizes(NAO_DJF_lt0_NAD_gt0) 
  print(NAO_DJF_lt0_NAD_gt0)


  slp_comp = dim_avg_n_Wrap(slp_DJF(NAD_NAO_DJF_gt5,:,:), 0)
  slp_comp2 = dim_avg_n_Wrap(slp_DJF(NAD_NAO_DJF_lt5,:,:), 0)
  slp_comp3 = dim_avg_n_Wrap(slp_DJF(NAO_DJF_gt5_NAD_lt0,:,:), 0)
  slp_comp4 = dim_avg_n_Wrap(slp_DJF(NAO_DJF_gt5_NAD_lt5,:,:), 0)
  slp_comp5 = dim_avg_n_Wrap(slp_DJF(NAO_DJF_gt0_NAD_lt0,:,:), 0)
  slp_comp6 = dim_avg_n_Wrap(slp_DJF(NAO_DJF_lt0_NAD_gt0,:,:), 0)
  ;*****************************t-test*******************************
  ;******************************************************************
  slp_std   = dim_variance_n_Wrap(slp_DJF(NAD_NAO_DJF_gt5,:,:), 0)
  slp_std   = sqrt(slp_std/nnumb)
  slp_std   = where(slp_std.eq.0,slp_std@_FillValue,slp_std)
  t_slp     = slp_comp/slp_std
  confi_slp = slp_comp
  confi_slp = student_t(t_slp, nnumb-1)

  slp_std2   = dim_variance_n_Wrap(slp_DJF(NAD_NAO_DJF_lt5,:,:), 0)
  slp_std2   = sqrt(slp_std2/nnumb2)
  slp_std2   = where(slp_std2.eq.0,slp_std2@_FillValue,slp_std2)
  t_slp2     = slp_comp2/slp_std2
  confi_slp2 = slp_comp2
  confi_slp2 = student_t(t_slp2, nnumb2-1)

  slp_std3   = dim_variance_n_Wrap(slp_DJF(NAO_DJF_gt5_NAD_lt0,:,:), 0)
  slp_std3   = sqrt(slp_std3/nnumb3)
  slp_std3   = where(slp_std3.eq.0,slp_std3@_FillValue,slp_std3)
  t_slp3     = slp_comp3/slp_std3
  confi_slp3 = slp_comp3
  confi_slp3 = student_t(t_slp3, nnumb3-1)

  slp_std4   = dim_variance_n_Wrap(slp_DJF(NAO_DJF_gt5_NAD_lt5,:,:), 0)
  slp_std4   = sqrt(slp_std4/nnumb4)
  slp_std4   = where(slp_std4.eq.0,slp_std4@_FillValue,slp_std4)
  t_slp4     = slp_comp4/slp_std4
  confi_slp4 = slp_comp4
  confi_slp4 = student_t(t_slp4, nnumb4-1)

  slp_std5   = dim_variance_n_Wrap(slp_DJF(NAO_DJF_gt0_NAD_lt0,:,:), 0)
  slp_std5   = sqrt(slp_std5/nnumb5)
  slp_std5   = where(slp_std5.eq.0,slp_std5@_FillValue,slp_std5)
  t_slp5     = slp_comp5/slp_std5
  confi_slp5 = slp_comp5
  confi_slp5 = student_t(t_slp5, nnumb5-1)

  slp_std6   = dim_variance_n_Wrap(slp_DJF(NAO_DJF_lt0_NAD_gt0,:,:), 0)
  slp_std6   = sqrt(slp_std6/nnumb6)
  slp_std6   = where(slp_std6.eq.0,slp_std6@_FillValue,slp_std6)
  t_slp6     = slp_comp6/slp_std6
  confi_slp6 = slp_comp6
  confi_slp6 = student_t(t_slp6, nnumb6-1)

;*********************************************************************
;*********************************************************************
  base=new(6,"graphic")
  plot=new(6,"graphic")

  Bres                =True
  Bres@gsnDraw        = False
  Bres@gsnFrame       = False 
  Bres@cnInfoLabelOn  = False
  Bres@gsnMaximize    = True
  ;Bres@gsnAddCyclic   = True
  Bres@gsnLeftString  =""
  Bres@gsnRightString =""

  Bres@mpFillOn       = False
  Bres@mpOutlineOn    = True
  Bres@mpGeophysicalLineThicknessF = "0.8"
  Bres@mpOutlineDrawOrder = "PostDraw"

  Bres@mpLimitMode    = "LatLon"
  Bres@mpMinLatF      = 0
  Bres@mpMaxLatF      = 80
  Bres@mpMinLonF      = -160
  Bres@mpMaxLonF      = 60
  Bres@mpCenterLonF   = 0
  Bres@mpShapeMode    = "FreeAspect"
  Bres@vpHeightF      = 0.5
  Bres@vpWidthF       = 0.8
  Bres@tmYLMode       = "Explicit"
  Bres@tmYLValues     = (/0,20,40,60,80/)
  Bres@tmYLLabels     = (/"EQ","20~S~o~N~N","40~S~o~N~N","60~S~o~N~N","80~S~o~N~N"/)
  Bres@tmXBMode="Explicit"
  Bres@tmXBValues=(/-150,-120,-90,-60,-30,0,30,60/)
  Bres@tmXBLabels=(/"150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0","30~S~o~N~E","60~S~o~N~E"/)
  
  Bres@gsnCenterString = "DJF-NAD&NAO.gt.0.5" ;子图的主标题
  Bres@gsnCenterStringFontHeightF   =0.035

  Bres@tmYROn        = False        ;关掉右边以及顶端的大刻度
  Bres@tmXTOn        = False

  Bres@tmXBMinorOn   = False   ;关掉所有的小刻度
  Bres@tmXTMinorOn   = False                  
  Bres@tmYLMinorOn   = False
  Bres@tmYRMinorOn   = False
 
  Bres@tmXBMajorLengthF = 0.008
  Bres@tmYLMajorLengthF = 0.008

  Bres@cnFillOn                 = True
  Bres@cnLinesOn                = True
  Bres@cnLineLabelsOn           = False
  Bres@lbLabelBarOn             = False
  Bres@cnLineLabelPlacementMode = "Constant"
  Bres@cnLineThicknessF         = 0.3
  Bres@cnLineDashPattern        = "0"      ;0为实线，1为虚线
  Bres@cnLineLabelBackgroundColor ="0"

 
  Bres2 = Bres
  Bres2@gsnCenterString = "DJF-NAD&NAO.lt.-0.5" ;子图的主标题

  Bres3 = Bres
  Bres3@gsnCenterString = "DJF-NAO.gt.0.5&NAD.lt.0" ;子图的主标题

  Bres4 = Bres
  Bres4@gsnCenterString = "DJF-NAO.gt.0.5&NAD.lt.0.5" ;子图的主标题

  Bres5 = Bres
  Bres5@gsnCenterString = "DJF-NAO.gt.0&NAD.lt.0" ;子图的主标题

  Bres6 = Bres
  Bres6@gsnCenterString = "DJF-NAO.lt.0&NAD.gt.0" ;子图的主标题

  base(0) = gsn_csm_contour_map_ce(wks, slp_comp, Bres)
  base(1) = gsn_csm_contour_map_ce(wks, slp_comp2, Bres2)
  base(2) = gsn_csm_contour_map_ce(wks, slp_comp3, Bres3)
  base(3) = gsn_csm_contour_map_ce(wks, slp_comp4, Bres4)
  base(4) = gsn_csm_contour_map_ce(wks, slp_comp5, Bres5)
  base(5) = gsn_csm_contour_map_ce(wks, slp_comp6, Bres6)

  sres                          = True                                                
  sres@gsnDraw                  = False                                                  
  sres@gsnFrame                 = False                                                
  sres@cnLineLabelsOn           = False                                                   
  sres@cnLinesOn                = False                                                   
  sres@cnInfoLabelOn            = False
  sres@cnFillOn                 = True   
  sres@gsnAddCyclic   			= True  
  sres@lbLabelBarOn             = False
  sres@cnInfoLabelOn            = False

  sres@cnFillDotSizeF           = 0.0005
  sres@cnLevelSelectionMode = "ExplicitLevels"  
  sres@cnLevels                 = 0.1
  sres@gsnAddCyclic             = True
  sres@cnFillPatterns           = 17

  plot(0) = gsn_csm_contour(wks, confi_slp, sres)
  plot(1) = gsn_csm_contour(wks, confi_slp2, sres)
  plot(2) = gsn_csm_contour(wks, confi_slp3, sres)
  plot(3) = gsn_csm_contour(wks, confi_slp4, sres)
  plot(4) = gsn_csm_contour(wks, confi_slp5, sres)
  plot(5) = gsn_csm_contour(wks, confi_slp6, sres)

  opt                           = True
  opt@gsnShadeFillType          = "pattern"
  opt@gsnShadeMid               = 17

  plot(0) = gsn_contour_shade(plot(0),-999.,0.1, opt)
  plot(1) = gsn_contour_shade(plot(1),-999.,0.1, opt)
  plot(2) = gsn_contour_shade(plot(2),-999.,0.1, opt)
  plot(3) = gsn_contour_shade(plot(3),-999.,0.1, opt)
  plot(4) = gsn_contour_shade(plot(4),-999.,0.1, opt)
  plot(5) = gsn_contour_shade(plot(5),-999.,0.1, opt)

  overlay(base(0), plot(0))
  overlay(base(1), plot(1))
  overlay(base(2), plot(2))
  overlay(base(3), plot(3))
  overlay(base(4), plot(4))
  overlay(base(5), plot(5))


  cnres                               = True
  cnres@gsnPanelLabelBar              = True
  cnres@gsnMaximize                   = True
  cnres@gsnPanelLabelPosition         = "Center"
  cnres@lbBoxEndCapStyle="TriangleBothEnds"   ;色标尾端样式
  cnres@lbLabelAlignment="BoxCenters"   ;色标标签对齐框中间

  gsn_panel(wks,base,(/3,2/),cnres)


  ; draw(base)
  ; frame(wks)


end
