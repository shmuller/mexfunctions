FC = gfortran

iiex_obj = iiex.o
ivex_obj = ivex.o
ixex_obj = ixex.o
xex1_obj = xex1.o bsplvb.o interv.o
xex2_obj = xex2.o bsplvb.o
xex3_obj = xex3.o bsplpp.o bsplvb.o ppvalu.o interv.o
xex4_obj = xex4.o bvalue.o interv.o
xiiex2_obj = xiiex2.o cubspl.o newnot.o
xiiex3_obj = xiiex3.o cubspl.o
xiiex4_obj = xiiex4.o bsplpp.o bsplvb.o
xiiiex1_obj = xiiiex1.o splint.o banfac.o banslv.o bsplpp.o bsplvb.o \
			  ppvalu.o interv.o round.o
xiiiex2_obj = xiiiex2.o splint.o banfac.o banslv.o bsplpp.o bsplvb.o \
			  newnot.o
xiiiex2m_obj = xiiiex2m.o splint.o banfac.o banslv.o bsplpp.o bsplvb.o \
			  newnot.o
xiiiex3_obj = xiiiex3.o titand.o splopt.o splint.o bvalue.o bsplvb.o \
			  banfac.o banslv.o interv.o
xivex1_obj = xivex1.o bsplpp.o bsplvb.o ppvalu.o interv.o smooth.o \
			  setupq.o chol1d.o
reinsch_obj = reinsch.o bsplpp.o bsplvb.o ppvalu.o interv.o smooth.o \
			  setupq.o chol1d.o
xivex2_obj = l2main.o setdatx2.o l2knts.o l2appr.o bsplvb.o \
			  bchfac.o bchslv.o bsplpp.o l2err.o ppvalu.o \
			  interv.o newnot.o
xivex3_obj = l2main.o setdatx3.o l2knts.o l2appr.o bsplvb.o \
			  bchfac.o bchslv.o bsplpp.o l2err.o ppvalu.o \
			  interv.o newnotfk.o
xivex4_obj = l2main.o setdatx4.o titand.o l2knts.o l2appr.o bsplvb.o \
			  bchfac.o bchslv.o bsplpp.o l2err.o ppvalu.o \
			  interv.o newnot.o
xvex_obj = xvex.o colloc.o colpnt.o difequ.o ppvalu.o interv.o knots.o \
			  eqblok.o putit.o bsplvd.o bsplvb.o bsplpp.o newnot.o \
			  slvblktext.o 
xviex2_obj = xviex2.o titand.o tautsp.o ppvalu.o interv.o
xviex3_obj = xviex3.o splint.o bsplvb.o banfac.o banslv.o \
			  bsplpp.o ppvalu.o interv.o
xviiex2_obj = xviiex2.o spli2d.o bsplvb.o banfac.o banslv.o interv.o \
			  bvalue.o 
xviiex3_obj = xviiex3.o spli2d.o bsplvb.o banfac.o banslv.o interv.o \
			  bspp2d.o ppvalu.o

targets = iiex ivex ixex xex1 xex2 xex3 xex4 xiiex2 xiiex3 xiiex4 \
          xiiiex1 xiiiex2 xiiiex2m xiiiex3 xivex1 \
	      xvex xviex2 xviex3 xviiex2 xviiex3

all : $(targets)

clean :
	rm $(targets) *.o

iiex : $(iiex_obj) 
	$(FC) $(iiex_obj) -o iiex

ivex : $(ivex_obj) 
	$(FC) $(ivex_obj) -o ivex

ixex : $(ixex_obj) 
	$(FC) $(ixex_obj) -o ixex

xex1 : $(xex1_obj) 
	$(FC) $(xex1_obj) -o xex1

xex2 : $(xex2_obj) 
	$(FC) $(xex2_obj) -o xex2

xex3 : $(xex3_obj)
	$(FC) $(xex3_obj) -o xex3

xex4 : $(xex4_obj)
	$(FC) $(xex4_obj) -o xex4

xiiex2 : $(xiiex2_obj)
	$(FC) $(xiiex2_obj) -o xiiex2

xiiex3 : $(xiiex3_obj)
	$(FC) $(xiiex3_obj) -o xiiex3

xiiex4 : $(xiiex4_obj)
	$(FC) $(xiiex4_obj) -o xiiex4

xiiiex1 : $(xiiiex1_obj) 
	$(FC) $(xiiiex1_obj) -o xiiiex1

xiiiex2 : $(xiiiex2_obj)
	$(FC) $(xiiiex2_obj) -o xiiiex2

xiiiex2m : $(xiiiex2m_obj)
	$(FC) $(xiiiex2m_obj) -o xiiiex2m

xiiiex3 : $(xiiiex3_obj) 
	$(FC) $(xiiiex3_obj) -o xiiiex3

xivex1 : $(xivex1_obj)
	$(FC) $(xivex1_obj) -o xivex1

reinsch : $(reinsch_obj)
	$(FC) $(reinsch_obj) -o reinsch

xivex2 : $(xivex2_obj)
	$(FC) $(xivex2_obj) -o xivex2

xivex3 : $(xivex3_obj)
	$(FC) $(xivex3_obj) -o xivex3

xivex4 : $(xivex4_obj)
	$(FC) $(xivex4_obj) -o xivex4

xvex : $(xvex_obj)
	$(FC) $(xvex_obj) -o xvex

xviex2 : $(xviex2_obj)
	$(FC) $(xviex2_obj) -o xviex2

xviex3 : $(xviex3_obj)
	$(FC) $(xviex3_obj) -o xviex3

xviiex2 : $(xviiex2_obj)
	$(FC) $(xviiex2_obj) -o xviiex2

xviiex3 : $(xviiex3_obj)
	$(FC) $(xviiex3_obj) -o xviiex3
