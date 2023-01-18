typedef struct {
  PIODOUBLE sig;
  PIODOUBLE listp[MAXSIMU];
  PIOFLOAT comap;
  PIOFLOAT dustmap;
  PIOFLOAT pdustmap;
  PIOFLOAT listofpix[MAXTHEOHPR];
  PIODOUBLE vspline[4];
  int istart;
  int iend;
  PIOINT ipix;
  PIOINT rg; //num ring
  PIOINT hrg;
  PIOFLOAT freefree;
  PIOINT adu;
  PIOFLOAT sadu;
  PIOFLOAT corr_nl;
  PIOFLOAT corr_cnn;
  PIOINT gi;
  PIODOUBLE dip;
  PIOFLOAT fsl; 
  PIOFLOAT phase;
  PIOFLOAT w;   //weights
  PIOFLOAT co;
  PIOFLOAT si;  // pixel pour 32 orientation et 12*32*32 pixel
  PIOBYTE  surv;
  PIOBYTE  ib;  // index bolo/detecteur
  PIOFLOAT hit; //hitcount
  PIOFLOAT wp;
  PIODOUBLE vi,vq,vu;
  PIODOUBLE model;
} hpix;

