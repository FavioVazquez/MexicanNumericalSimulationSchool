LDFLAGS = -O2    -fconvert=little-endian   
FFLAGS  = -O2    -fconvert=little-endian   
FC = gfortran

BIN     = PMstartM PMmodels PMpowerM PMmainM PM_to_ASCIIm
OBJS    = PMauxiliary.o PMstartM.o PMpowerM.o PMmainM.o\
            PM_to_ASCIIm.o PM_to_ASCIIm.o 
HEADERS = PMparameters.h PMinitial.h PMlists.h 


# Para el caso de un modelo LCDM, se deben compilar estos programas

PMmainM: PMmainM.o PMauxiliary.o
	$(FC) $(LDFLAGS) -o PMmainM PMmainM.o PMauxiliary.o $(LIB)

PMpowerM: PMpowerM.o PMauxiliary.o
	$(FC) $(LDFLAGS) -o PMpowerM PMpowerM.o PMauxiliary.o $(LIB)


PMmodels: PMmodels.o PMauxiliary.o
	$(FC) $(LDFLAGS) -o PMmodels PMmodels.o $(LIB)

PMstartM: PMstartM.o PMparameters.h
	$(FC) $(LDFLAGS) -o  PMstartM PMstartM.o  $(LIB)

# Estos son algunos programas auxiliares
#PMstartM.o: PMstartM.f PMparameters.h

PM_to_ASCIIm: PM_to_ASCIIm.o PMauxiliary.o
	$(FC) $(LDFLAGS) -o   PM_to_ASCIIm PM_to_ASCIIm.o PMauxiliary.o  $(LIB)

#PMauxiliary.o: PMparameters.h 

.f90.o: 
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f $(BIN) $(OBJS) 
