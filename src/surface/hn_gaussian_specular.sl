//Gaussian function of mean 0 and width beta
float gaussian (float beta,x)
{
        return((1/(beta*sqrt(2*PI)))*exp(-x*x/(2*beta*beta)));
}

surface hn_gaussian_specular(float roughness = 0.3; color specular_color = 1;)
{
	vector Nf = faceforward(normalize(N),I);	
	vector R = reflect(I,Nf);
	vector Rn = normalize(R);
	color spec = 0;
	float theta_r = acos(Rn.Nf);
	
	illuminance(P,Nf,PI/2)
	{
		//Calculate incident angle of incoming light
		vector Ln = normalize(L);	
		float theta_i = acos(Ln.Nf);
		spec += Cl*gaussian(roughness,theta_i-theta_r);
	}
	Ci = specular_color*spec*Os;
	Oi = Os;
}		
