displacement d_RP_simplenoise_disp(
	float Km = 0.03;
	float noiFreq = 10;)
{
	/*Local variables*/
	float noi;

	/*Compute a noise field on the shader coordinate sys*/
	noi = noise(transform("shader",P*noiFreq));

	/*Modify the position of P and recalculate the normal */
	P += normalize(N)*noi*Km;
	N = calculatenormal(P);
} 
