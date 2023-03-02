import numpy as np
import os, sys
import matplotlib.pyplot as plt
import healpy as hp
import getopt

#=================================================================================
# INITIALIZE FoCUS class
#=================================================================================
import foscat.Synthesis as synthe


def usage():
    print(' This software is a demo of the foscat library:')
    print('>python demo.py -n=8 [-c|--cov][-s|--steps=3000][-x|--xstat][-p|--p00][-g|--gauss][-k|--k5x5]')
    print('-n : is the nside of the input map (nside max = 256 with the default map)')
    print('--cov (optional): use scat_cov instead of scat.')
    print('--steps (optional): number of iteration, if not specified 1000.')
    print('--xstat (optional): work with cross statistics.')
    print('--p00   (optional): Loss only computed on p00.')
    print('--gauss (optional): convert Venus map in gaussian field.')
    print('--k5x5  (optional): Work with a 5x5 kernel instead of a 3x3.')
    exit(0)
    
def run(signal,nside,rank):
  
    if rank == 0 : 
      print('[Foscat]')
  
  
    cov=False
    #nside=-1
    nstep=100
    docross=True
    dop00=False
    dogauss=False
    KERNELSZ=3   
    
    print(len(signal))


       

    if cov:
        import foscat.scat_cov as sc
        print('Work with ScatCov')
    else:
        import foscat.scat as sc
        print('Work with Scat')
    
    np.save("map_test.npy",signal)
    
    return 0;
  
    imap = signal
    imap[imap==hp.UNSEEN] = 0;
    #=================================================================================
    # Function to reduce the data used in the FoCUS algorithm 
    #=================================================================================
    def dodown(a,nside):
        nin=int(np.sqrt(a.shape[0]//12))
        if nin==nside:
            return(a)
        return(np.mean(a.reshape(12*nside*nside,(nin//nside)**2),1))
      
    #=================================================================================
    # DEFINE A PATH FOR scratch data
    # The data are storred using a default nside to minimize the needed storage
    #=================================================================================
    scratch_path = '/export/home/tfoulquier/FOSCAT/scratch'
    outname='TEST_model'

    proj ='I'

    idx=hp.nest2ring(nside,np.arange(12*nside**2)) 
    im = hp.read_map('/export/home1/tfoulquier/out_maps/test_rstep2_353psb_nside%s_%s.fits'%(nside,proj),nest=True)[idx]*1000
   
    #imap = hp.read_map('/export/home1/tfoulquier/out_maps/test_rstep2_857GHz_ODUST_GRAD_nside%s_%s.fits'%(nside,proj),nest=True)[idx]

    #=================================================================================
    # Get data
    #=================================================================================
    
    if dogauss:
        idx=hp.ring2nest(nside,np.arange(12*nside*nside))
        idx1=hp.nest2ring(nside,np.arange(12*nside*nside))
        cl=hp.anafast(im[idx])
        im=hp.synfast(cl,nside)[idx1]

        hp.mollview(im,cmap='jet',nest=True)
        plt.show()
    

    lam=1.2
    if KERNELSZ==5:
        lam=1.0
    #=================================================================================
    # COMPUTE THE WAVELET TRANSFORM OF THE REFERENCE MAP
    #=================================================================================
    scat_op=sc.funct(NORIENT=4,          # define the number of wavelet orientation
                     KERNELSZ=KERNELSZ,  # define the kernel size
                     OSTEP=-1,           # get very large scale (nside=1)
                     LAMBDA=lam,
                     TEMPLATE_PATH=scratch_path,
                     slope=1.0,
                     gpupos=0,
                     all_type='float32')
      
    #=================================================================================
    # DEFINE A LOSS FUNCTION AND THE SYNTHESIS
    #=================================================================================
    #F(353,353) = F(x,x)
    def lossX(x,scat_operator,args):
        
        ref = args[0]
        im  = args[1]


        if docross:
            learn=scat_operator.eval(im,image2=x,Imaginary=True)
        else:
            learn=scat_operator.eval(x)
            
        if dop00:
            loss=scat_operator.bk_reduce_mean(scat_operator.bk_square(ref.P00[0,0,:]-learn.P00[0,0,:]))
        else:
            loss=scat_operator.reduce_sum(scat_operator.reduce_mean(scat_operator.square(ref-learn)))
            

        return(loss)
    
    
    
    if docross:
        refX=scat_op.eval(im,image2=im,Imaginary=True)
    else:
        refX=scat_op.eval(im)
    

    loss1=synthe.Loss(lossX,scat_op,refX,im)
    sy = synthe.Synthesis([loss1])

    #=================================================================================
    # RUN ON SYNTHESIS
    #================================================================================= 
    omap=sy.run(imap,
                DECAY_RATE=0.999,
                NUM_EPOCHS = nstep,
                LEARNING_RATE = 0.3,
                EPSILON = 1E-7)

    #=================================================================================
    # STORE RESULTS
    #=================================================================================
    if docross:
        start=scat_op.eval(im,image2=imap)
        out =scat_op.eval(im,image2=omap)
    else:
        start=scat_op.eval(imap)
        out =scat_op.eval(omap)

    nitt = 1
    np.save('out_foscat_test/out_model_test_%s_map_%d.npy'%(proj,nside),omap)    

    np.save('out_foscat_test/in_model_test_%s_map_%d_itt%d.npy'%(proj,nside,nitt),im)
    np.save('out_foscat_test/st_model_test_%s_map_%d_itt%d.npy'%(proj,nside,nitt),imap)
    np.save('out_foscat_test/out_model_test_%s_map_%d_itt%d.npy'%(proj,nside,nitt),omap)
    np.save('out_foscat_test/out_model_test_%s_log_%d_itt%d.npy'%(proj,nside,nitt),sy.get_history())

    refX.save('out_foscat_test/in_model_test_%s_%d_itt%d'%(proj,nside,nitt))
    start.save('out_foscat_test/st_model_test_%s_%d_itt%d'%(proj,nside,nitt))
    out.save('out_foscat_test/out_model_test_%s_%d_itt%d'%(proj,nside,nitt))

    print('Computation Done')
    sys.stdout.flush()

#im= np.load('map_test.npy')
#run(im,'plop',64,0)

