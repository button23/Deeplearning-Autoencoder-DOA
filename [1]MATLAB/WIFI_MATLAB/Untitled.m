
a = wlanConstellationMap(interleaver_data,Nbpscs);
b = wlanConstellationDemap(a,0,Nbpscs,'soft');
c= abs(a);

% Constellation Demapper(algorithm)
llrAggregate = zeros(Ncbps, Nsym);
A = sqrt(1/170);
qamBits= zeros(Nbpscs,52);
for oo = 1  :  Nsym
        for ll = 1: length(a) 
            rr = real(a(ll,oo));
            im = imag(a(ll,oo));
            llr(1) = rr; 
            llr(2) = - abs(rr) + 8 * A;
            llr(3) = - abs( abs(rr) - 8 * A ) + 4 * A;
            llr(4) = - abs( abs( abs(rr) - 8 * A ) - 4 * A) + 2 * A;
            
            llr(5) = im; 
            llr(6) = - abs(im) + 8 * A;
            llr(7) = - abs( abs(im) - 8 * A ) + 4 * A;
            llr(8) = - abs( abs( abs(im) - 8 * A ) - 4 * A) + 2 * A;
            
            qamBits(:,ll) = llr * (abs(H_VHT(ll)))^2 ; 
%             qamBits(:,ll) = llr ; 


        end
        llrAggregate(:,oo) = reshape(qamBits,Ncbps,1);        
end

demod_Algorithm = llrAggregate > 0;
numError = biterr(demod_Algorithm,interleaver_data)
ber = numError / (Ncbps*Nsym)

compp = reshape(qamBits,416,1);
b_scale = b/10^8;

ratio = ( b_scale ./ compp);

