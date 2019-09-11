load('simulation_28_2.mat');
met_name='NAA';
B1=eval([met_name,'_err',num2str(1)]);
SNR=SNR./3.4;
figure();
for kk=1:5
    for i=1:5
        A=eval([met_name,'_val',num2str(i)]);
        if i==1
            B=zeros(size(A));
        else
            B=eval([met_name,'_err',num2str(i-1)]);
        end
        if i==1
            minA=min(min(A))- 2*max(max(B1));
            maxA=max(max(A))+ 2*max(max(B1));
        end
        subplot(5,5,5*(kk-1)+i);
        errorbar(SNR,A(kk,:),B(kk,:),'-o','MarkerSize',2,'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',2);
        axis([0 110./3.5 minA maxA]);  
        id=round([100 65 50 30 20 10]./3.45) +1;
        set(gca, 'XTick',fliplr(id), 'FontSize',1)
        set(gca,'XDir','reverse');
        if (kk==1)
            if i==1
                title(['accl=',num2str(accl_factor(i)),'(CSI)']);
            else
                title(['accl=',num2str(accl_factor(i))]);
            end
        end
        if (kk==3 && i==1)
            ylabel([met_name,' peak area']);
        end
        if (kk==5 && i==3)
            xlabel(['SNR']);
        end
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',9,'FontWeight','Bold', 'LineWidth', 2);
    end
end
dim = [.03 .05 .1 .1];
str = 'B0 map 5';
annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
dim = [.03 .23 .1 .1];
str = 'B0 map 4';
annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
dim = [.03 .4 .1 .1];
str = 'B0 map 3';
annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
dim = [.03 .58 .1 .1];
str = 'B0 map 2';
annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
dim = [.03 .75 .1 .1];
str = 'B0 map 1';
annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
% figure();
% for kk=1:5
%     for i=1:1
%         A=eval([met_name,'_val',num2str(i)]);
%         B=eval([met_name,'_err',num2str(i)]);
%         if i==1
%             minA=min(min(A))- 2*max(max(B));
%             maxA=max(max(A))+ 2*max(max(B));
%         end
%         subplot(5,1,1*(kk-1)+i);
%         errorbar(SNR,A(kk,:),'-o','MarkerSize',2,'MarkerEdgeColor','red','MarkerFaceColor','red','linewidth',2);
% %         axis([0 100 minA maxA]);  
%         if (kk==1)
%             title(['CSI']);
%         end
%         if (kk==3 && i==1)
%             ylabel([met_name,' peak area']);
%         end
%         if (kk==5 && i==1)
%             xlabel(['SNR']);
%         end
%         set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',9,'FontWeight','Bold', 'LineWidth', 2);
%     end
% end
% dim = [.03 .05 .1 .1];
% str = 'B0 map 1';
% annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
% dim = [.03 .23 .1 .1];
% str = 'B0 map 2';
% annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
% dim = [.03 .4 .1 .1];
% str = 'B0 map 3';
% annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
% dim = [.03 .58 .1 .1];
% str = 'B0 map 4';
% annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
% dim = [.03 .75 .1 .1];
% str = 'B0 map 5';
% annotation('textbox',dim,'String',str,'FontSize',10,'FontWeight','Bold');
close all;
met_name='NAA';
metabolite=cell(0);
metabolite{1}='NAA';metabolite{2}='Lac';metabolite{3}='Cr';metabolite{4}='Cho';metabolite{5}='Myo';
figure();
for kk=1:1
    for i=1:5
            B1=eval([metabolite{kk},'_err',num2str(1)]);
            A=eval([metabolite{kk},'_val',num2str(i)]);
            if i==1
                B=zeros(size(A));
            else
                B=eval([metabolite{kk},'_err',num2str(i-1)]);
            end
            if i==1
                minA=min(min(A))- 2*max(max(B1));
                maxA=max(max(A))+ 2*max(max(B1));
            end
            if i>0
            subplot_tight(1,5,5*(kk-1)+(i));
            ThreeDBarWithErrorBars(A,B,B1,minA+.1,maxA);
            zoom(1)
            set(gca,'FontName','Arial','FontSize',8,'FontWeight','Bold',  'LineWidth', 1.5);
    %         set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',9,'FontWeight','Bold', 'LineWidth', 2);
            if kk==1 
                if i==1
                    title(['accl=',num2str(accl_factor(i)),'(CSI)']);
                else
                    title(['accl=',num2str(accl_factor(i)),'(RAPID)']);
                end
            end
            if i==0
                xlabel('SNR');
                ylabel(['\DeltaB_{' int2str(0) '} map']);
                zlabel([metabolite{kk},' area']);
            end
                set(gca,'XTickLabel',([30 20 12 8 6 4.5 2.5]));
                set(gca,'YTick',([1:1:5]));
                set(gca,'YTickLabel',(['1'; '2'; '3'; '4'; '5' ]));
            end
    end
end
