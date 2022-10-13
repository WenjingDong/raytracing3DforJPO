function [] = savedata(fnum, nfile, t, kxvec, flag)
     
     disp(['write data ' num2str(fnum) ' of ' num2str(nfile)]);
     if flag==0
         fileID = fopen('output.bin','w');
     else
         fileID = fopen('output.bin','a');
     end
     fwrite(fileID, kxvec, 'double');
     fclose(fileID);
     
     if flag==0
         fileID = fopen('time.bin','w');
     else
         fileID = fopen('time.bin','a');
     end
     fwrite(fileID, t, 'double');
     fclose(fileID);
end