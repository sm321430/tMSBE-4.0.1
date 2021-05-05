for j=1:49
for k=1:2
x(j,k)=CAVsim(j).scaledInt(k);
end
end


theta_int = @(theta,ext_ind, cav_ind) asin(sin(theta)*ext_ind/cav_ind);
