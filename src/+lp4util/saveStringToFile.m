function saveStringToFile(s, file, mode)

if nargin < 3
    mode = 'wt'; % write as text
end

fid = fopen(file, mode);
fprintf(fid, '%s', s);
fclose(fid);

end

