//Gaussian function of mean 0 and width beta
float gaussian (float beta,x)
{
        return((1/(beta*sqrt(2*PI)))*exp(-x*x/(2*beta*beta)));
}

surface marschner_fur(float beta_R = 0.3, alpha_R = 0.1745, eta = 1.55;
		color specular_color = 1;)
{
	//Constants
	const float beta_TT = beta_R/2;
	const float beta_TRT = 2*beta_R;
	const float alpha_TT = -alpha_R/2;
	const float alpha_TRT = -3*alpha_R/2;
	
	//Forward Normal
	vector Nf = faceforward(normalize(N),I);
	
	//basis vectors for shading coordinate system
	vector T = normalize(dPdv); //tangent vector along hair direction
	vector W = normalize(Nf^T); //second basis vector in normal plane
	vector V = -normalize(I);

	//Normal plane projections of incident,scatterd vectors 
	vector L_phi,V_phi;
	V_phi = normalize((V.W)*W+(V.Nf)*Nf);
	
	//angles in shading coordinate system
	float Phi_i,Theta_i,Phi_r,Theta_r,Theta_h;
	Phi_r = acos(V_phi.W);
	Theta_r = acos(V_phi.V);
	
	//BRDF components
	float M_R = 0.0;
	
	//Total scattered value
	color S = 0.0;

	illuminance(P,Nf,PI/2)
	{
		//Normalized light direction vector
		vector Ln = normalize(N);	
		L_phi = normalize((Ln.W)*W + (Ln.Nf)*Nf);
		//Normalized light direction vector
		Phi_i = acos(L_phi.W);
		Theta_i = acos(Ln.L_phi);
		Theta_h = (Theta_i + Theta_r)/2;
		//TODO: gamma_i,sigma_a,h
		
		float Phi = Phi_r - Phi_i;
		//Bravais' indices
		float eta_1 = sqrt(eta^2 - (sin(Theta_i))^2 ) / cos(Theta_i);
		float eta_2 = eta * cos(Theta_i) / sqrt(1 - eta^-2* (sin(Theta_i))^2 );
		
		S += ((M_R(beta_R,alpha_R,Theta_h) * N_R(eta_1,eta_2,getGamma_i(Phi,eta_1,0),getDphidh(0,eta_1,getGamma_i(Phi,eta_1,0))))+
		
			(M_TT(beta_TT,alphaTT,Theta_h) * N_TT(eta_1,eta_2,getGamma_i(Phi,eta_1,1),sigma_a,h,getDphidh(1,eta_1,h)))+
			
			(M_TRT(beta_TRT,alpha_TRT,Theta_h) * N_TRT(eta_1,eta_2,getGamma_i(Phi,eta_1,2),gamma_t,sigma_a,h)))
			*Cl;
	}
	Ci = S*Os;
	Oi = Os;	
}

float getDphidh(int p, float eta_1,float gamma_i)
{
	float h = getH(gamma_i);
	return ( (2*p/eta_1) * (1.0/sqrt(1.0-(h^2/eta_1^2))) - (2/sqrt(1.0 - h^2)) );
}

float getGamma_i(float Phi,eta_1; int p)
{
	//TODO: solve cubic equation 5 to get gamma_i
}

float N_R(float eta_1,eta_2,gamma_i,dPhidh)
{
	return fresnel(eta_1,eta_2,gamma_i) / (2* dPhidh(0,eta_1,getH(gamma_i)))
}

float getH(float gamma_i)
{
	return arcsin(gamma_i);
}

float N_TRT(float eta_1,eta_2, gamma_i,gamma_t,sigma_a,h)
{
	//TODO
}

float M_TRT(float beta_TRT,alpha_TRT,Theta_h)
{
	//TODO
}

float N_TT(float eta_1,eta_2,gamma_i,sigma_a,h,dPhidh)
{
	//TODO
}

float M_TT(float beta_TT,alpha_TT,Theta_h)
{
	//TODO
}


float getDphidh()
{
	//TODO
}
float M_R(float beta_R, float alpha_R, float Theta_h)
{
	return gaussian(beta_R,Theta_h-alpha_R)
}	
