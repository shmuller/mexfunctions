#INC=-I$(MDSPLUS_DIR)/include
#LIB=-L$(MDSPLUS_DIR)/lib -lMdsIpShr

#EXE=test_mdsconnect

INC=-I"C:\PROGRA~1\MDSPLUS\DEVTOOLS\include"
LIB="C:\PROGRA~1\MDSPLUS\DEVTOOLS\lib\mdsipshr.lib" -lws2_32

EXE=test_mdsconnect.exe

$(EXE):	test_mdsconnect.c
	gcc test_mdsconnect.c -o $(EXE) $(INC) $(LIB)
