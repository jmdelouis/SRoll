pro donewtf,rk,nrk,ibstart=ibstart,dowrite=dowrite,doapp=doapp

if not keyword_set(ibstart) then ibstart=0

bol=['353-3a','353-3b','353-5a','353-5b','353-4a','353-4b','353-6a','353-6b','353-1','353-2','353-7','353-8']
calib=[3.48625e-14,3.37349e-14,2.86904e-14,2.83174e-14,3.27368e-14,3.23329e-14,2.34751e-14,2.12574e-14,5.70346e-14,6.20394e-14,4.70167e-14,4.42892e-14]
bol2=['23_353_3a','24_353_3b','53_353_5a','54_353_5b','32_353_4a','33_353_4b','63_353_6a','64_353_6b','05_353_1','13_353_2','45_353_7','85_353_8']

gg=[[    0.0045420751,    0.0029670298,    0.0014776677], $
    [    0.0040155819,    0.0026851201,    0.0015305540], $
    [    0.0021077225,    0.0017270296,    0.0015807816], $
    [    0.0051315009,    0.0041933192,    0.0023808829], $
    [    0.0021543254,    0.0025556159,    0.0021544109], $
    [    0.0046122191,    0.0031805323,    0.0015515991], $
    [    0.0057770136,    0.0040070780,    0.0025883108], $
    [    0.0061261651,    0.0048697178,    0.0037163580], $
    [    0.0031208132,    0.0027873084,    0.0019119671], $
    [    0.0022460001,    0.0024911353,    0.0017339265], $
    [    0.0042414795,    0.0027782916,    0.0014256445], $
    [    0.0038798249,    0.0031748097,    0.0018297135]]

bol=['353-2','353-7','353-8']
calib=[6.20394e-14,4.70167e-14,4.42892e-14]
bol2=['13_353_2','45_353_7','85_353_8']

gg=[[    0.0022460001,    0.0024911353,    0.0017339265], $
    [    0.0042414795,    0.0027782916,    0.0014256445], $
    [    0.0038798249,    0.0031748097,    0.0018297135]]

;gg2=reform(read_binary('/redtruck/delouis/M4RESULT/RD12_RC2/353_SWB/VECT/353-3a_offsets_PROD_REP6_ns_RD12_353GHz_X2',data_typ=5),8,41)
;gg=reform(read_binary('/redtruck/delouis/M4RESULT/RD12_RC2/353_SWB/VECT/353-3a_offsets_PROD_REP6_RD12_353GHz_X2',data_typ=5),12,41)
;gg(0:7,*)=gg2
;print,pioreadringindexgrp(b,e,n,pioopentoigrp('/data/dmc/MISS03/DATA/calTOIs','r'))
b=read_binary('/redtruck/delouis/PROD_RD12_RC4/DATA/BEGINRINGINDEX',data_type=14)
e=read_binary('/redtruck/delouis/PROD_RD12_RC4/DATA/ENDRINGINDEX',data_type=14)

nstep=262144LL
RINGSIZE=27664LL
frq=262144/180./60. 
beg=lonarr(5)
beg(0)=0.5*frq
dstep=1
for i=1,4 do begin
    beg[i]=beg[i-1]+frq*dstep   
    dstep*=2
end

xref=(beg(1:*)+beg(0:*))/2
xx=shift(dindgen(nstep)-nstep/2,nstep/2)

maxpix=785000LL

na=['0','1','2','3']

for ib=ibstart,ibstart+nrk-1 do begin ;n_elements(bol)-1 do begin
    scale=1.0
    bolo=bol[ib]
    bolo2=bol2[ib]
    a=reform(read_binary('/redtruck/delouis/M4RESULT/RD12_RC2/353_SWB/VECT/353-3a_offsets_PROD_REP6_RD12_353GHz_X2',data_typ=5),12,41)

    xx=(2^(dindgen(3)+1)+2^(dindgen(3)+2)-1)/2

    res=dblarr(100,100)
    amp=dindgen(100)*3E-2
    tau=dindgen(100)*2E-1
    xxx=shift(dindgen(nstep)-nstep/2,-nstep/2)
    for i=0,99 do for j=0,99 do res(i,j)=total((a(ib,32:34)-float(amp[i]/(complex(tau[j]*frq,xx*frq,/double))))^2)
    ;!P.multi=[0,1,2]
    ;contour,-res,nlev=600,/fill
    err=min(res,p)
    aa=amp[p mod 100]
    bb=tau[p/100]
    ;plot,xx*frq,a(ib,32:34),xr=[0,256],yr=[0,0.01]
    filt1=float(1/(complex(4*frq,xxx,/double)))
    filt2=float(1/(complex(8*frq,xxx,/double)))
    ;oplot,xxx,aa*filt1,col=250
    ;oplot,xxx,aa*filt2,col=80
    print,aa,bb
    ;filt(*)=complex(0,0)
    ;for i=1,50 do filt(frq*(i-0.5):frq*(i+0.5)-1)=complex(1,0)*float(aa/(complex(bb*frq,i*frq,/double)))
    ;filt(beg[2]:beg[3]-1)=complex(1,0)
    ;filt(nstep-beg[3]+1:nstep-beg[2])=complex(1,0)

    if keyword_set(dowrite) then begin
        if keyword_set(doapp) then begin
            openu,u1,'/redtruck/delouis/PROD_RD12_RC4/PBR_JMD/'+bolo+'_REP6_TT1',/get_lun,/append
            openu,u2,'/redtruck/delouis/PROD_RD12_RC4/PBR_JMD/'+bolo+'_REP6_TT2',/get_lun,/append
        end else begin
            openw,u1,'/redtruck/delouis/PROD_RD12_RC4/PBR_JMD/'+bolo+'_REP6_TT1',/get_lun
            openw,u2,'/redtruck/delouis/PROD_RD12_RC4/PBR_JMD/'+bolo+'_REP6_TT2',/get_lun
        end
    end

    if keyword_set(dowrite) then begin
        if not keyword_set(doapp) then begin
            vide=fltarr(RINGSIZE)
            for ir=0,239 do begin
                writeu,u1,vide
                writeu,u2,vide
            end
        end
    end
    ;print,piocreatepoiobject('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_TT1','PIOFLOAT')
    ;print,piocreatepoiobject('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_TT2','PIOFLOAT')

    for ir=rk,26050 do begin
        com='begin='+string(b[ir]-nstep/4)+';end='+string(b[ir]+nstep/2*8+nstep/2)
        
        ;idx=pioread('/data/dmc/MISS03/DATA/e2e_common_TOI/'+bolo+'_HPRIDX_ABER_TotalFlag_dx11',com=com)
        ;toi=pioread('/data/dmc/MISS03/DATA/calTOIs/'+bolo2+'_LFER4_JC_v61',com=com)
        ;rg=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_hit',com='ring='+string(ir))
                                ;ph=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_ptg',com='ring='+string(ir))
                                ;th=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_ptg_TUPLE_1',com='ring='+string(ir))
        idx=readdmc('/data/dmc/MISS03/DATA/e2e_common_TOI/'+bolo+'_HPRIDX_ABER_TotalFlag_dx11',b[ir]-nstep/4,b[ir]+nstep/2*8+nstep/2)
        toi=readdmc('/data/dmc/MISS03/DATA/calTOIs/'+bolo2+'_LFER4_JC_v61',b[ir]-nstep/4,b[ir]+nstep/2*8+nstep/2)
        rg=readdmc('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_hit',ir*RINGSIZE,(ir+1)*RINGSIZE-1)
        ph=readdmc('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_ptg',ir*RINGSIZE,(ir+1)*RINGSIZE-1)
        th=readdmc('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_ptg_TUPLE_1',ir*RINGSIZE,(ir+1)*RINGSIZE-1)

        toi-=median(toi)
        tois=fltarr(n_elements(toi))
        tois2=fltarr(n_elements(toi))
        for ii=0,7 do begin
            x1=ii*nstep/2
            x2=ii*nstep/2+nstep-1
            tf=fft(toi(x1:x2),-1)
            tf2=filt1*tf
            res=fft(tf2,1)
            tois((x1+x2)/2-nstep/4:(x1+x2)/2+nstep/4)=res(nstep/2-nstep/4:nstep/2+nstep/4)
            tf2=filt2*tf
            res=fft(tf2,1)
            tois2((x1+x2)/2-nstep/4:(x1+x2)/2+nstep/4)=res(nstep/2-nstep/4:nstep/2+nstep/4)
        end
        res=dblarr(RINGSIZE)
        res2=dblarr(RINGSIZE)
        hres=dblarr(RINGSIZE)
        endidx=e[ir]-b[ir]-1
        if endidx ge maxpix then  endidx=maxpix
        for i=nstep/4,endidx+nstep/4 do if idx(i) ge 0 then hres(idx(i))+=1
        for i=nstep/4,endidx+nstep/4 do if idx(i) ge 0 then res(idx(i))+=tois(i)
        for i=nstep/4,endidx+nstep/4 do if idx(i) ge 0 then res2(idx(i))+=tois2(i)
        
        if total(hres gt 0) then begin
            res(where(hres gt 0))/=hres(where(hres gt 0))*calib[ib]
            res2(where(hres gt 0))/=hres(where(hres gt 0))*calib[ib]

            ;!p.multi=[0,1,2]
            ;plot,res*scale
            ;plot,res2*scale
            
            ;stop
            if keyword_set(dowrite) then begin
                
                writeu,u1,float(res*scale)
                writeu,u2,float(res2*scale)
            end

            ;err=piowritepoiobject(res*scale,'/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_TT1','PIOFLOAT','ring='+string(ir))
            ;err=piowritepoiobject(res2*scale,'/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_TT2','PIOFLOAT','ring='+string(ir))
           
            print,bolo,ir,err,total((hres-rg)^2)
            
            nRG=27664
        
            ;rr=fltarr(nRG,5)
            ;rr(*,0)=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_R0',com='ring='+string(ir))
            ;rr(*,1)=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_R1',com='ring='+string(ir))
            ;rr(*,2)=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_R2',com='ring='+string(ir))
            ;rr(*,3)=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_R3',com='ring='+string(ir))
            ;rr(*,4)=pioread('/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_R4',com='ring='+string(ir))
            
            
            ;mat=dblarr(5,5)
            ;t=where(hres gt 0,nn)
            ;for i=0,4 do for j=0,4 do mat(i,j)=total(rr(t,i)*rr(t,j))
            
            ;vec=fltarr(5)
            ;for j=0,4 do vec(j)=total(scale*res(t)*rr(t,j))
            
            ;cc=cramer(mat,vec)
            ;print,cc
;            if scale eq 1 then begin
;                scale=sqrt(total(a(ib,32:34)*cc(1:3))/total(cc(1:3)*cc(1:3)))
;                err=piowritepoiobject(res*scale,'/data/dmc/MISS03/DATA/PBR_JMD/'+bolo+'_REP6_TF2','PIOFLOAT','ring='+string(ir))
;            end
            ;oplot,xx*frq,cc(1:*)*scale,col=80
            ;stop
        end

        if total((hres-rg)^2) gt 1 then begin
            print,'PBS RING ',bolo,ir,total((hres-rg)^2)
            stop
            return
        end
    end
    
    if keyword_set(dowrite) then begin
        free_lun,u1
        free_lun,u2
    end
    
end


end
