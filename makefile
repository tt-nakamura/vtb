NTL = -lntl -lgmp
OBJ = ConvolCode.o RSCode.o HCCode.o GF2Xlib.o

main: main.o $(OBJ)
	g++ main.o $(OBJ) $(NTL)

fig4: fig4.o $(OBJ)
	g++ fig4.o $(OBJ) $(NTL)