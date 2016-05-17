//Gaussian function of mean 0 and width beta
float gaussian (float beta,x)
{
        return((1/(beta*sqrt(2*PI)))*exp(-x*x/(2*beta*beta)));
}

float N_R(float Phi)
{
	return cos(Phi/2);
}

float N_TT(float gamma_TT,Phi)
{
	return gaussian(gamma_TT,PI-Phi);
}

float N_TRT(float Phi, gamma_g,I_g)
{
	extern point P;
	float G_angle = (float cellnoise(transform("world",P)))* 0.26179 + 0.52358;
	return ( cos(Phi/2) +I_g*gaussian(gamma_g,G_angle - Phi));
}

float M_TRT(float beta_TRT,alpha_TRT,Theta_h)
{
	return gaussian(beta_TRT,Theta_h-alpha_TRT);
}

float M_TT(float beta_TT,alpha_TT,Theta_h)
{
	return gaussian(beta_TT,Theta_h-alpha_TT);
}

float M_R(float beta_R, alpha_R, Theta_h)
{
	return gaussian(beta_R,Theta_h-alpha_R);
}	

surface hair_sadeghi(float beta_R = 0.3, alpha_R = 0.1745, gamma_TT=0.3, gamma_g=0.3, I_g=0.5,I_R=1.0,I_TT=1.0,I_TRT=1.0,I_diff=0.5;
		color specular_color = 1; output varying color R = 0, TT = 0, TRT = 0, diffuse = 0;)
{
	//Dummy variables
	color _R = 0;
	color _TT = 0;
	color _TRT = 0;
	color _diff = 0;
	
	//Constants
	uniform float beta_TT = beta_R/2.0;
	uniform float beta_TRT = 2*beta_R;
	uniform float alpha_TT = -alpha_R/2;
	uniform float alpha_TRT = -3*alpha_R/2;
	
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
		vector Ln = normalize(L);	
		L_phi = normalize((Ln.W)*W + (Ln.Nf)*Nf);

		Phi_i = acos(L_phi.W);
		Theta_i = acos(Ln.L_phi);
		Theta_h = (Theta_i + Theta_r)/2.0;
		float Phi = Phi_r - Phi_i;
		
		_diff += I_diff*(Ln.Nf)*Cl*specular_color*Os;
		_R   += I_R*M_R(beta_R,alpha_R,Theta_h) * N_R(Phi)*Cl*specular_color*Os;
		_TT  += I_TT*M_TT(beta_TT,alpha_TT,Theta_h) * N_TT(gamma_TT,Phi)*Cl*specular_color*Os;
		_TRT += I_TRT*M_TRT(beta_TRT,alpha_TRT,Theta_h) * N_TRT(Phi,gamma_g,I_g)*Cl*specular_color*Os;
		S += (I_R*(M_R(beta_R,alpha_R,Theta_h) * N_R(Phi))+
		
			(I_TT*M_TT(beta_TT,alpha_TT,Theta_h) * N_TT(gamma_TT,Phi))+
			
			(I_TRT*M_TRT(beta_TRT,alpha_TRT,Theta_h) * N_TRT(Phi,gamma_g,I_g)))
			*Cl;
	}
	R = _R;
	TT = _TT;
	TRT = _TRT;
	diffuse = _diff;
	Ci = _R+_TT+_TRT+_diff;
	Oi = Os;	
}