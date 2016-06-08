function KL = kl_div_gaussians(mu0,sig0,mu1,sig1)

KL = 0.5*(trace(sig1\sig0)+(mu1-mu0)'*sig1\(mu1-mu0)-size(sig0,1)+log(det(sig1))-log(det(sig0)));