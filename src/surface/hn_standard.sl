surface hn_standard(color specular_color=1, diffuse_color=1; float specular_factor=1.0, roughness=0.3, diffuse_factor=1)
{
	normal Nf = normalize(faceforward(N,I));
	vector V = -normalize(I);
	color result = 0;
	color spec,diff = 0;

	illuminance(P,Nf,PI/2)
	{
		vector Ln = normalize(L);
		vector H = normalize(Ln+V);
		float hdotn = H.Nf;	
		spec = Cl*pow(max(hdotn,0),1/roughness);
		diff = Cl*(Ln.Nf);
		result += spec*specular_factor*specular_color
				 +diff*diffuse_factor*diffuse_color; 
	}
	Ci = Os*result;
	Oi = Os;
}	
		
