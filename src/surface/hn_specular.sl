surface hn_specular(float specular_factor=1.0, roughness=0.3)
{
	normal Nf = normalize(faceforward(N,I));
	vector V = -normalize(I);
	color spec = 0;

	illuminance(P,Nf,PI/2)
	{
		vector Ln = normalize(L);
		vector H = normalize(Ln+V);
		float hdotn = H.Nf;	
		spec += Cl*pow(max(hdotn,0),1/roughness); 
	}
	Ci = Os*spec*specular_factor;
	Oi = Os;
}	
		
