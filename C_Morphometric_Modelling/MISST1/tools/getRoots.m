function protocol = getRoots(protocol)
% function protocol = getRoots(protocol)
% adds the bessel function roots necessary for diffusion signal computation
% to the protocol structure
%

protocol.roots_cyl = BesselJ_RootsCyl(40);
protocol.roots_sphere = BesselJ_RootsSphere(40);
protocol.roots_plane = (1:length(protocol.roots_cyl))*2-1;
end
