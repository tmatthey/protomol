%XYZ Reader
%
%[x] xyzread(filename);
function [x,x1,x2,x3,x4,x5] = xyzread(filename)

     % read header
     n=0;
     newline=sprintf('\n');
     fid = fopen(filename);
     while(n<2)
       c=fread(fid,1,'char');
       if(c==newline)
         n=n+1;
       end
     end

     % read positions
     n = 1;
     x=[];
     x1=[];
     x2=[];
     x3=[];
     x4=[];
     x5=[];
     map =[];
     map = setstr(map);

     fprintf('XYZ Reader ');
     while(isempty(ferror(fid)))
       clear name;
       name=fscanf(fid, '%s', [1,1]);
       p=fscanf(fid, '%g %g %g', [1,3]);
       if(~isempty(p))
          a=[name '          '];
          a=a(1:10);
          i=1;
          while(i<=size(map,1))
            if(strcmp(map(i,:),a))
               break;
            end;
            i=i+1;
          end
          if(i > size(map,1))
               map(size(map,1)+1,:)=a;
               i = size(map,1);
          end
          x(n,1:3)=p;
          x(n,4)=i;
          if(i==1)
            x1(size(x1,1)+1,1:3) = x(n,1:3);
          elseif (i==2)
            x2(size(x2,1)+1,1:3) = x(n,1:3);
	  elseif (i==3) 
            x3(size(x3,1)+1,1:3) = x(n,1:3);
	  elseif (i==4) 
            x4(size(x4,1)+1,1:3) = x(n,1:3);
	  else 
            x5(size(x5,1)+1,1:3) = x(n,1:3);
          end

%          fprintf('%d:%s %g %g %g %d \n',n,name,x(n,1),x(n,2),x(n,3),x(n,4));
          if(mod(n,1000)==0)
            fprintf('.');
          end
          n=n+1;
       end
     end;

     fclose(fid);

     if(n>1000)
       fprintf(' ');
     end
     fprintf('read %d position(s) and %d type(s) (',n-1,size(map,1));
     for i=1:size(map,1)
       for j=1:size(map,2)
         if(map(i,j)~=' ')
           fprintf('%s',map(i,j));
         end
       end
       if(i<size(map,1))
         fprintf(',');
       end
     end
     fprintf(').\n');
