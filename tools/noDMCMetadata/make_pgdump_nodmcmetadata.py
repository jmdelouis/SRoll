#!/usr/bin/env python

import os
import sys
import glob
import stools
import time


DUMPS = ["/redtruck/SimuData/PGDUMP_20171213_Wednesday", # M3
         "/pscratch1/RD12_data/PGDUMP_20171213_Wednesday", # M4
         "/sps/planck/SimuData/RD12_data/PGDUMP_20171213_Wednesday", # CC
         "/scratch/cnt0028/ias1717/SHARED/RD12_data/PGDUMP_20171213_Wednesday"] # occigen


def usage_and_abort():
  print "usage: %s [--verbose] [--nowrite] [--replace] <full_object_name>" % sys.argv[0]
  print "       <object(s)>: full object name(s), accepts wildcards"
  print "       [--verbose]: if set, will display more info to ease debuging"
  print "       [--nowrite]: if set, will not try to create the no_dmc_metadata.txt file"
  print "       [--replace]: if set, will overwrite existing no_dmc_metadata.txt file"
  exit(1)


# COPY dmc_users (user_id, user_name, user_mail, user_groups, user_defgroup, firstname, lastname) FROM stdin;
DMCUSERS="""
10000	dmcadmin	support.l2@planck.fr	{1}	1
101	jaatroko	jha@kurp.hut.fi	{6}	6	Juha	Aatrokoski
102	ekeihane	elina.keihanen@helsinki.fi	{6}	6	Elina	Keihanen
103	rkeskita	rejio.keskitalo@gmail.com	{3,2}	3	Reijo	Keskitalo
104	sgratton	stg20@cam.ac.uk	{3,2}	3	Steven	Gratton
106	sosborne	sjo32@stanford.edu	{3,2}	3	Stephen	Osborne
107	jstarck	jstarck@cea.fr	{3,2}	3	Jean-Luc	Starck
108	gsavini	gs@star.ucl.ac.uk	{3,2}	3	Giorgio	Savini
109	jaumont	jonathan.aumont@ias.u-psud.fr	{3,2}	3	Jonathan	Aumont
10	pdasilva	dasilva@iap.fr	{3,2}	3	Pierre	Da Silva
110	cnorth	chris.north@astro.cf.ac.uk	{3,2}	3	Christopher	North
112	ygiraudh	Yannick.Giraud-Heraud@apc.univ-paris-diderot.fr	{3,2}	3	Yannick	Giraud-Heraud
113	cchiang	hsinc@Princeton.edu	{3,2}	3	Cynthia	Chiang
114	abanday	Anthony.Banday@irap.omp.eu	{3,2}	3	Anthony	Banday
115	lveltz	lionel.veltz@ias.u-psud.fr	{3,2}	3	Lionel	Veltz
116	naghanim	nabila.aghanim@ias.u-psud.fr	{3,2}	3	Nabila	Aghanim
117	odor	olivier.p.dore@jpl.nasa.gov	{3,2}	3	Olivier	Dore
118	jborrill	jdborrill@lbl.gov	{3,2}	3	Julian	Borrill
119	kgorski	Krzysztof.M.Gorski@jpl.nasa.gov	{3,2}	3	Krzysztof	Gorski
120	smitra	sanjit@iucaa.ernet.in	{3,2}	3	Sanjit	Mitra
121	dnovikov	d.novikov@imperial.ac.uk	{3,2}	3	Dmitri	Novikov
122	jdelabro	delabrouille@apc.univ-paris7.fr	{3,2}	3	Jacques	Delabrouille
123	hpeiris	h.peiris@ucl.ac.uk	{3,2}	3	Hiranya	Peiris
124	jmelin	jean-baptiste.melin@cea.fr	{3,2}	3	Jean-Baptiste	Melin
125	pdeberna	paolo.debernardis@roma1.infn.it	{3,2}	3	Paolo	de Bernardis
126	smasi	silvia.masi@roma1.infn.it	{3,2,4}	3	Silvia	Masi
129	jgrain	julien.grain@ias.u-psud.fr	{3,2}	3	Julien	Grain
12	cmercier	claude.mercier@ias.u-psud.fr	{3,2}	3	Claude	Mercier
132	mkunz	Martin.Kunz@unige.ch	{3,2}	3	Martin	Kunz
133	gluzzi	gluzzi@lal.in2p3.fr	{3,2}	3	Gemma	Luzzi
134	shildebr	srh@caltech.edu	{3,2}	3	Sergi	Hildebrandt
135	bvantent	Bartjan.Van-Tent@th.u-psud.fr	{3,2}	3	Bartjan	van Tent
137	rkneissl	rkneissl@alma.cl	{3,2}	3	Ruediger	Kneissl
138	alewis	antony@cosmologist.info	{3,2}	3	Antony	Lewis
139	tjageman	Thomas.Jagemann@rssd.esa.int	{3,2}	3	Thomas	Jagemann
13	coxborro	oxborrow@space.dtu.dk	{3,2}	3	Carol Anne	Oxborrow
140	lspencer	locke.spencer@astro.cf.ac.uk	{3,2}	3	Locke	Spencer
141	mbridges	mb435@mrao.cam.ac.uk	{3,2}	3	Michael	Bridges
142	abenoitl	benoitl@iap.fr	{3,2}	3	Aurelien	Benoit-Levy
143	ylongval	yuying.longval@ias.u-psud.fr	{3,2}	3	Yuying	Longval
145	hdole	herve.dole@ias.u-psud.fr	{3,2}	3	Herve	Dole
146	cahernan	chm@cefca.es	{3,2}	3	Carlos	Hernandez-Monteagudo
147	dsutton	sutton@ast.cam.ac.uk	{3,2}	3	David	Sutton
148	dhanson	duncan.hanson@gmail.com	{3,2}	3	Duncan	Hanson
149	jdunkley	j.dunkley@physics.ox.ac.uk	{3,2}	3	Joanna	Dunkley
14	crenault	rcecile@in2p3.fr	{3,2,4}	3	Cecile	Renault
150	pcarvalh	fpmnpc2@cam.ac.uk	{3,2}	3	Pedro	Carvalho
151	ghurier	hurier@lpsc.in2p3.fr	{3,2}	3	Guillaume	Hurier
152	pbielewi	Pawel.Bielewicz@sissa.it	{3,2}	3	Pawel	Bielewicz
153	tdechele	dechelet@iap.fr	{3,2}	3	Typhaine	Dechelette
154	fnoviell	fabio.noviello@manchester.ac.uk	{3,2}	3	Fabio	Noviello
155	mremazei	Remazeil@apc.univ-paris-diderot.fr	{3,2}	3	Mathieu	Remazeilles
156	jlesgour	lesgourg@lapp.in2p3.fr	{3,2}	3	Julien	Lesgourgues
157	fpaci	fpaci@apc.univ-paris7.fr	{3,2}	3	Francesco	Paci
158	mbrown01	mbrown@jb.man.ac.uk	{3,2}	3	Michael	Brown
159	aminiu01	antoine.miniussi@ias.u-psud.fr	{3,2}	3	Antoine	Miniussi
15	crosset	cyrille.rosset@apc.univ-paris-diderot.fr	{3,2}	3	Cyrille	Rosset
160	gcastex	Guillaume.Castex@apc.univ-paris-diderot.fr	{3,2}	3	Guillaume	Castex
161	mtucci	tucci@lal.in2p3.fr	{3,2}	3	Marco	Tucci
162	lmaurin	loic.maurin@apc.univ-paris-diderot.fr	{3,2}	3	Loic	Maurin
163	cgauthie	cgauthie@apc.univ-paris7.fr	{3,2}	3	Christopher	Gauthier
164	carmitag	armitage-caplan@physics.ox.ac.uk	{3,2}	3	Charmaine	Armitage-Caplan
165	mfrommer	mona.frommert@unige.ch	{3,2}	3	Mona	Frommert
167	lsanselm	sanselme@lpsc.in2p3.fr	{3,2}	3	Lilian	Sanselme
168	efalgaro	edith.falgarone@lra.ens.fr	{3,2}	3	Edith	Falgarone
169	afraisse	afraisse@princeton.edu	{3,2}	3	Aurelien	Fraisse
16	dclement	d.clements@imperial.ac.uk	{3,2}	3	David	Clements
170	stechene	techene@iap.fr	{3,2}	3	Sibylle	Techene
171	dmortloc	d.mortlock@imperial.ac.uk	{3,2}	3	Daniel	Mortlock
172	flacasa	fabien.lacasa@ias.u-psud.fr	{3,2}	3	Fabien	Lacasa
173	silic	stephane.ilic@ias.u-psud.fr	{3,2}	3	Stephane	Ilic
174	jgudmund	jgudmund@princeton.edu	{3,2}	3	Jon	Gudmundsson
176	tjaffe	tess.jaffe@irap.omp.eu	{3,2}	3	Tess	Jaffe
177	sdrouot	sebastien.drouot@lpsc.in2p3.fr	{6}	6	Sebastien	Drouot
178	lpagano	luca.pagano@roma1.infn.it	{3,2}	3	Luca	Pagano
179	felsner	elsner@iap.fr	{3,2}	3	Franz	Elsner
17	dfleisch	dsf@ast.cam.ac.uk	{3,2}	3	Dominique	Fleischmann
180	vguillet	vincent.guillet@ias.u-psud.fr	{3,2}	3	Vincent	Guillet
181	fsureau	florent.sureau@cea.fr	{3,2}	3	Florent	Sureau
182	cparisel	camille.parisel@cesr.fr	{3,2}	3	Camille	Parisel
183	khuffenb	huffenbe@physics.miami.edu	{3,2}	3	Kevin	Huffenberger
184	fboulang	francois.boulanger@ias.u-psud.fr	{3,2}	3	Francois	Boulanger
185	lcolom01	colombo@usc.edu	{3,2}	3	Loris	Colombo
186	lknox	lknox@physics.ucdavis.edu	{3,2}	3	Lloyd	Knox
187	mmillea	mmillea@ucdavis.edu	{3,2}	3	Marius	Millea
188	zhou	hou@ucdavis.edu	{3,2}	3	Zhen	Hou
189	hsangher	hss@ast.cam.ac.uk	{3,2}	3	Hardip	Sanghera
190	jbowyer	j.bowyer07@imperial.ac.uk	{3,2}	3	Jude	Bowyer
191	arahlin	arahlin@princeton.edu	{3,2}	3	Alexandra	Rahlin
192	malves	marta.alves@ias.u-psud.fr	{3,2}	3	Marta	Alves
193	mlanger	mathieu.langer@ias.u-psud.fr	{3,2}	3	Mathieu	Langer
194	lfauvet	lfauvet@rssd.esa.int	{3,2}	3	Lauranne	Fauvet
195	tghosh	tuhin.ghosh@ias.u-psud.fr	{3,2}	3	Tuhin	Ghosh
196	bwandelt	benwandelt@googlemail.com	{3,2}	3	Benjamin D.	Wandelt
197	ccombet	celine.combet@lpsc.in2p3.fr	{3,2}	3	Celine	Combet
198	pserra	paolo.serra@ias.u-psud.fr	{3,2}	3	Paolo	Serra
199	bruizgra	bearg@ugr.es	{3,2}	3	Beatriz	Ruiz-Granados
19	dharriso	dlh@ast.cam.ac.uk	{3,2}	3	Diana L	Harrison
1	acatalan	catalano@lpsc.in2p3.fr	{3,2,4}	3	Andrea	Catalano
20000	opsman	opsman@planck.fr	{2,3}	2	Planck Operations	Manager
200	fnati	federico.nati@roma1.infn.it	{3,2}	3	Federico	Nati
201	bracine	benjamin.racine@apc.univ-paris7.fr	{3,2}	3	Benjamin	Racine
202	forieux	orieux@iap.fr	{3,2}	3	Francois	Orieux
203	ghuey	ggh1@cosmology.name	{3,2}	3	Greg	Huey
204	psutter	sutter@iap.fr	{3,2}	3	Paul	Sutter
205	bcomis	comis@lpsc.in2p3.fr	{3,2}	3	Barbara	Comis
206	acollier	acollier@lbl.gov	{3,2}	3	Aaron	Collier
207	hjimenez	jimenez@iap.fr	{3,2}	3	Hugo	Jimenez-Perez
208	lvignon	lvignon@ias.u-psud.fr	{3,2}	3	Ludivine	Vignon
209	amangill	anna.mangilli@gmail.com	{3,2}	3	Anna	Mangilli
20	eglez	eglez@ast.cam.ac.uk	{3,2}	3	Eduardo	Gonzalez Solares
210	ksmith01	kmsmith@astro.princeton.edu	{3,2}	3	Kendrick	Smith
211	flevrier	francois.levrier@ens.fr	{3,2}	3	Francois	Levrier
213	gefstath	gpe@ast.cam.ac.uk	{3,2}	3	George	Efstathiou
214	mroman	mroman@apc.univ-paris7.fr	{3,2}	3	Matthieu	Roman
215	dalina	Dana.Alina@irap.omp.eu	{3,2}	3	Dana	Alina
216	dpietrob	davide.pietrobon@jpl.nasa.gov	{3,2}	3	Davide	Pietrobon
217	radam	adam@lpsc.in2p3.fr	{3,2}	3	Remi	Adam
218	ecalabre	Erminia.Calabrese@astro.ox.ac.uk	{3,2}	3	Erminia	Calabrese
219	ganiano	gonzalo.aniano@ias.u-psud.fr	{3,2}	3	Gonzalo	Aniano
21	ehivon	hivon@iap.fr	{3,2}	3	Eric	Hivon
221	bbertinc	benjamin.bertincourt@ias.u-psud.fr	{3,2}	3	BENJAMIN	BERTINCOURT
222	mspinell	spinelli@lal.in2p3.fr	{3,2}	3	Marta	Spinelli
223	sgalli	gallis@iap.fr	{3,2}	3	Silvia	Galli
22	epointec	etienne.pointecouteau@irap.omp.eu	{3,2}	3	Etienne	Pointecouteau
24	fbouchet	bouchet@iap.fr	{3,2,4}	3	Francois R.	Bouchet
25	fcouchot	couchot@lal.in2p3.fr	{3,2}	3	Francois	Couchot
26	fdesert	Francois-Xavier.Desert@obs.ujf-grenoble.fr	{3,2,4}	3	Francois-Xavier	Desert
27	fferoz	f.feroz@mrao.cam.ac.uk	{6}	6	Farhan	Feroz
29	fpajot	francois.pajot@ias.u-psud.fr	{3,2,4}	3	Francois	Pajot
2	achambal	antoine.chamballu@cea.fr	{3,2}	3	Antoine	Chamballu
30000	hfiproduct	opsman@planck.fr	{5}	5
30	fpiacent	Francesco.Piacentini@roma1.infn.it	{3,2,4}	3	Francesco	Piacentini
31	ftouze	touze@lal.in2p3.fr	{3,2}	3	Francois	Touze
32	glagache	guilaine.lagache@ias.u-psud.fr	{3,2}	3	Guilaine	Lagache
34	gpatanch	patanchon@apc.univ-paris7.fr	{3,2}	3	Guillaume	Patanchon
35	gprezeau	Gary.Prezeau@jpl.nasa.gov	{3,2}	3	Gary	Prezeau
36	grocha	graca@caltech.edu	{3,2}	3	Graca	Rocha
37	jbartlet	James.Bartlett@apc.univ-paris-diderot.fr	{3,2}	3	James	Bartlett
38	jbernard	jean-philippe.bernard@cesr.fr	{3,2}	3	Jean-Philippe	Bernard
39	jcardoso	cardoso@tsi.enst.fr	{3,2}	3	Jean-Francois	Cardoso
3	acoulais	alain.coulais@obspm.fr	{3,2}	3	Alain	Coulais
40000	lfiman	opsman@planck.fr	{2,3,4}	4
40	jcolley	colley@apc.univ-paris7.fr	{3,2}	3	Jean-Marc	Colley
41	jdelouis	delouis@iap.fr	{3,2}	3	Jean-Marc	Delouis
42	jlamarre	jean-michel.lamarre@obspm.fr	{3,2,4}	3	Jean-Michel	Lamarre
43	jmacias	macias@lpsc.in2p3.fr	{3,2,4}	3	Juan Francisco	Macias-Perez
44	jpuget	jean-loup.puget@ias.u-psud.fr	{3,2,4}	3	Jean-Loup	Puget
45	jtauber	jtauber@rssd.esa.int	{3,2}	3	Jan	Tauber
46	kbenabed	benabed@iap.fr	{3,2}	3	Karim	Benabed
47	kdassas	Karin.Dassas@ias.u-psud.fr	{3,2}	3	Karin	Dassas
48	kganga	ganga@apc.univ-paris7.fr	{3,2}	3	Ken	Ganga
4	ajaffe	a.jaffe@imperial.ac.uk	{3,2}	3	Andrew	Jaffe
50	lmontier	ludovic.montier@cesr.fr	{3,2,4}	3	Ludovic	Montier
51	lperotto	perotto@lpsc.in2p3.fr	{3,2}	3	Laurence	Perotto
52	lvibert	laurent.vibert@ias.u-psud.fr	{3,2}	3	Laurent	Vibert
53	mashdown	maja1@mrao.cam.ac.uk	{3,2}	3	Mark	Ashdown
54	mbucher	bucher@apc.univ-paris7.fr	{3,2}	3	Martin	Bucher
55	mdouspis	marian.douspis@ias.u-psud.fr	{3,2}	3	Marian	Douspis
56	mgiard	martin.giard@irap.omp.eu	{3,2,4}	3	Martin	Giard
57	lejeune	lejeune@apc.univ-paris7.fr	{3,2}	3	Maude	Le Jeune
58	mpiat	piat@apc.univ-paris7.fr	{3,2,4}	3	Michel	Piat
59	mreineck	martin@mpa-garching.mpg.de	{3,2}	3	Martin	Reinecke
5	amoneti	moneti@iap.fr	{3,2}	3	Andrea	Moneti
60	mtristra	tristram@lal.in2p3.fr	{3,2}	3	Matthieu	Tristram
61	nponthie	nicolas.ponthieu@ias.u-psud.fr	{3,2}	3	Nicolas	Ponthieu
62	oherent	herent@iap.fr	{3,2}	3	Olivier	Herent
63	operdere	perderos@lal.in2p3.fr	{3,2,4}	3	Olivier	Perdereau
65	priant	riant@iap.fr	{3,2}	3	Philippe	Riant
66	rstompor	radek@apc.univ-paris7.fr	{3,2}	3	Radek	Stompor
67	scolombi	colombi@iap.fr	{3,2}	3	Stephane	Colombi
68	smottet	mottet@iap.fr	{3,2}	3	Sylvain	Mottet
69	splaszcz	plaszczy@lal.in2p3.fr	{3,2}	3	Stephane	Plaszczynski
6	bcrill	Brendan.P.Crill@jpl.nasa.gov	{3,2}	3	Brendan	Crill
70	sprunet	prunet@iap.fr	{3,2}	3	Simon	Prunet
71	sroubero	rouberol@iap.fr	{3,2}	3	Stephane	Rouberol
72	sschurr	schurr@ipac.caltech.edu	{3,2}	3	Steve	Schurr
73	tmatsumu	tm@astro.caltech.edu	{6}	6	Tomotake	Matsumura
74	vstolyar	vlad@ast.cam.ac.uk	{3,2}	3	Vladislav	Stolyarov
76	aducout	ducout@iap.fr	{3,2}	3	Anne	Ducout
79	wjones	wcjones@princeton.edu	{3,2}	3	William	Jones
7	brusholm	rusholme@ipac.caltech.edu	{3,2}	3	Ben	Rusholme
81	fboudol	boudol@apc.univ-paris7.fr	{6}	6	Florian	Boudol
82	rgael	roudier@apc.univ-paris7.fr	{3,2}	3	Gael	Roudier
84	tpears01	tjp@astro.caltech.edu	{3,2}	3	Tim	Pearson
85	thandley	thh@ipac.caltech.edu	{3,2}	3	Tom	Handley
87	asauve	alexandre.sauve@cesr.fr	{3,2}	3	Alexandre	Sauve
88	fmelot	frederic.melot@lpsc.in2p3.fr	{3,2}	3	Frederic	Melot
89	ebrelle	breelle@apc.univ-paris7.fr	{3,2}	3	Eric	Breelle
90	oforni	olivier.forni@irap.omp.eu	{3,2}	3	Olivier	Forni
91	pcamus	camus@grenoble.cnrs.fr	{3,2}	3	Philippe	Camus
92	gguyot	guy.guyot@ias.u-psud.fr	{3,2}	3	Guy	Guyot
93	gpoullea	gilles.poulleau@ias.u-psud.fr	{3,2}	3	Gilles	Poulleau
94	vhervier	Veronique.Hervier@ias.u-psud.fr	{3,2}	3	Veronique	Hervier
95	mmiville	marc-antoine.miville-deschenes@ias.u-psud.fr	{3,2}	3	Marc-Antoine	Miville-Deschenes
96	mcharra	maryse.charra@ias.u-psud.fr	{3,2}	3	Maryse	Charra
97	tmaciasz	thierry.maciaszek@oamp.fr	{3,2}	3	Thierry	Maciaszek
98	wlee	wplee@ipac.caltech.edu	{3,2}	3	Wen-Piao	Lee
99	jsygnet	sygnet@iap.fr	{3,2}	3	Jean-Francois	Sygnet
9	cfilliar	filliarc@lal.in2p3.fr	{3,2}	3	Clement	Filliard
"""

user_dict = dict()
for user in DMCUSERS.split("\n"):
  if user.strip() == "":
    continue
  user_fields = user.split("\t")
  user_dict[user_fields[0]] = user_fields[1]

if len(sys.argv) == 1:
  usage_and_abort()

PGDUMP = None
for tmp in DUMPS:
  if os.path.exists( tmp):
    PGDUMP = tmp
    print "Using %s\n" % PGDUMP
    break
if PGDUMP == None:
  print "PostgreSQL dump file not found: "
  print DUMPS

params = sys.argv[1:]

optVerbose = False
optNowrite = False
optReplace = False
for i, objname in enumerate( params):
  if objname.strip().lower() == "--verbose":
    del params[i]
    optVerbose = True
  elif objname.strip().lower() == "--nowrite":
    del params[i]
    optNowrite = True
  elif objname.strip().lower() == "--replace":
    del params[i]
    optReplace = True
  elif objname.startswith("-"):
    print "error: unkown option " + objname
    print
    usage_and_abort()
assert not (optNowrite and optReplace)

for param in params:
  for objname in glob.glob( param):
    metadataname = objname + "/no_dmc_metadata.txt"
    print metadataname
    if os.path.exists( metadataname) and (not optReplace):
      print "already exists, skipping\n"
      continue
    nodmctime = time.strftime("%Y-%m-%d %H:%M:%S")
    dump = stools.subprocess_check_output( 'grep -P "%s\t" %s' % (objname, PGDUMP))
    lines = dump.split("\n")
    if len( lines) < 2:
      print "grep PGDUMP: not enough entries found: "
      print dump
      exit( 0)

    for line in lines:
      cols = line.strip().split("\t")
      if optVerbose:
        print cols
      if len( cols) == 11:
#         Table "public.dmc_backendobject"
#      Column       |       Type        | Modifiers
#-------------------+-------------------+-----------
# object_id         | bigint            | not null
# parent_id         | bigint            | not null
# name              | character varying | not null
# backendid         | bigint            | not null
# backendmetaname   | character varying | not null
# backendname       | character varying | not null
# backendalias      | character varying | not null
# backendtype       | character varying | not null
# nativeobject      | character varying | not null
# transfertfunction | bigint            | not null
# tobebackup        | character varying | not null
        object_id = cols[0]
        backendid = cols[3]
        backendmetaname = cols[4]
        backendname = cols[5]
        backendtype = cols[7]
        assert backendtype.endswith("Object")
        backendtype = backendtype[:-len("Object")]
        assert backendmetaname == objname

      if ((backendtype == "POI") and (len( cols) == 16)) or len( cols) == 15:
#      Table "public.dmc_piolib_toiobject"
#     Column     |       Type        | Modifiers
#----------------+-------------------+-----------
# object_id      | bigint            | not null
# parent_id      | bigint            | not null
# datatype       | character varying | not null
# imoref         | bigint            | not null
# physicalproc   | character varying | not null
# memoryresident | character varying | not null
# softdeleted    | character varying | not null
# ioversion      | integer           | not null
# hostname       | character varying | not null
# backendname    | character varying | not null
# swapendian     | integer           | not null
# iooffset       | bigint            | not null
# flagchunksize  | bigint            | not null
# filesize       | bigint            | not null
# bc_hexa        | character varying | not null
        datatype = cols[2]
        iooffset = cols[11]
        flagchunksize = cols[12]
        assert int( flagchunksize) == 16384
        assert int( iooffset) == 128

    dump = stools.subprocess_check_output( 'grep $%s %s' % (backendid, PGDUMP))
    lines = dump.split("\n")
    if len( lines) < 2:
      print "grep PGDUMP: not enough entries found: "
      print dump
      exit( 0)

    for line in lines:
      cols = line.strip().split("\t")
      if optVerbose:
        print cols
      if len( cols) == 5:
#   Table "public.dmc_piolib_index"
#     Column      |  Type  | Modifiers
#-----------------+--------+-----------
# object_id       | bigint | not null
# parent_id       | bigint | not null
# backendobjectid | bigint | not null
# beginindex      | bigint | not null
# endindex        | bigint | not null
        begidx = cols[3]
        endidx = int(cols[4]) - 1

      if len( cols) == 10:
#                                            Table "public.dmc_objects"
#      Column       |            Type             |                            Modifiers
#-------------------+-----------------------------+-----------------------------------------------------------------
# object_id         | bigint                      | not null default nextval('dmc_objects_object_id_seq'::regclass)
# parent_id         | bigint                      | not null default nextval('dmc_objects_parent_id_seq'::regclass)
# type_id           | bigint                      | not null
# author            | bigint                      | not null
# group_id          | bigint                      | not null
# mask              | character varying           | not null
# creation_date     | timestamp without time zone | not null
# modification_date | timestamp without time zone | not null
# version           | bigint                      | not null
# label             | character varying           |
        author = cols[3]
        creadate = cols[6]
        modifdate = cols[7]

        if author in user_dict:
          author = user_dict[author]
        else:
          author = author + " (make_pgdump_nodmcmetada.py)"

    metadatatxt = """noDMCmetadata_date : {nodmctime}
noDMCmetadata_version : 0.0.2
TOItype : {backendtype}
Datatype : {datatype}
BeginIndex : {begidx}
EndIndex : {endidx}
Author : {author}
Date : {creadate}
Keyword list: 0 keyword(s) found
Location : unknown
Backendname : {backendname}
Flagchunksize : {flagchunksize}
Iooffset : {iooffset}
""".format(**locals())
    print metadatatxt

    if not optNowrite:
      if not os.path.exists( objname):
        print "%s directory not found, no_dmc_metadata.txt not written" % objname
      else:
        try:
          with open( metadataname, "w") as f:
            f.write( metadatatxt)
          print "written to disk\n"
        except:
          print "Error when writing to %s, check your rights" % metadataname

