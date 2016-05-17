surface hn_diffuse(float Ka=0.1, Kd=0.5, wrap=1., gam=1)
{
	normal Nf = normalize(faceforward(N,I));
	color Cwd=0;
	float wrp = 1-0.5*wrap;
	float wa = acos(wrp*2-1); //wrap angle
	
	illuminance(P, Nf, wa)
	{
		vector Ln = normalize(L);
		float diffuze;
		float dotp = 0.5*(1+Ln.Nf);
		if (dotp<-wrp)
		{
			diffuze = pow((wrp-dotp)/(wrp),gam);
		}
		else
			diffuze = pow((dotp-wrp)/(1-wrp),gam);
		Cwd += Cl*diffuze;
	}
	Ci = Cs*Os*(Ka*ambient() + Kd*Cwd);
	Oi = Os;
}
