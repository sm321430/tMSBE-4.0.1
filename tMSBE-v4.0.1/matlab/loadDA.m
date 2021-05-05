function a = loadDA(name,size)
% Read an 2d array of doubles 
% Name is name of file
% Size is [time-steps, K-steps]

if numel(size)~=2
   disp('loadDA(): Cannot read more or less than 2d arrays')
   name
   size
end

fid = fopen(name);
a = fread(fid,size(2:-1:1),'double');
fclose(fid);

a = a';

end
