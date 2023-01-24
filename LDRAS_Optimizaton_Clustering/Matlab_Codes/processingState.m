%% FUNCTION TO KNOW THE STATE OF THE PROCESS
function it_current=processingState(l,j,it_current)

    if l==1
        disp('0 %');
        it_current = 0;
    else
        iter = l/j;
        switch(it_current)
            case 0
                if(iter >= 0.1)
                    it_current = 10;
                    disp('10 %');
                end      
            case 10
                if(iter >= 0.2)
                    it_current = 20;
                    disp('20 %');
                end                
            case 20
                if(iter >= 0.3)
                    it_current = 30;
                    disp('30 %');
                end                
            case 30
                if(iter >= 0.4)
                    it_current = 40;
                    disp('40 %');
                end               
            case 40
                if(iter >= 0.5)
                    it_current = 50;
                    disp('50 %');
                end                    
            case 50
                if(iter >= 0.6)
                    it_current = 60;
                    disp('60 %');
                end                   
            case 60
                if(iter >= 0.7)
                    it_current = 70;
                    disp('70 %');
                end                       
            case 70
                if(iter >= 0.8)
                    it_current = 80;
                    disp('80 %');
                end                     
            case 80
                if(iter >= 0.9)
                    it_current = 90;
                    disp('90 %');
                end                       
            case 90
                if(iter >= 1.0)
                    it_current = 100;
                    disp('100 %');
                end
        end
    end
end