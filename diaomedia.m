function B = diaomedia(BOLD_timeries)
j=size(BOLD_timeries,2);   %360 脑区数  时间点x脑区数
B=cell(j,1);

for m=1:j%j
    M=BOLD_timeries(:,m);        
   for x=1:j
    X=BOLD_timeries(:,x);
    for y=1:j
        Y=BOLD_timeries(:,y);
        if (m==x || m==y || x==y)
            m_value = 0;
        else
        [bp,p] = xiugai(X, Y, M); 
        if ((p(1,1)<0.05 && p(1,2)<0.05) || p(1,5)<0.05)
                if (p(1,3)>=0.05) 
                    m_value = 1;
                end
         if (p(1,3)<0.05) 
            if (bp(1,5) * bp(1,3) > 0)
                m_value = bp(1,5) / bp(1,4);  %ab/c
            else
                m_value = 0;
            end
         end
        end                
        B{m,1}(x,y)=m_value;
        end
    end
  end 
end
end

