function[]=plot_index_data_standardized(X_raw,num_fig,GMModel)
%!! X_raw must be standardized
figure(num_fig)
clf

[~,sy]=size(X_raw);

for i=1:sy
    for j=1:sy
        if i<=j
            subplot(sy,sy,(i-1)*sy+j)
            if i==j
                figure(num_fig)
                hold on
                [nb_clust,~]=size(GMModel.mu);
                sum_pdf=zeros(1,length(min(X_raw(:,i)):0.1:max(X_raw(:,i))));
                for k=1:nb_clust
                    my_mu=GMModel.mu(k,i);
                    my_sig=sqrt(GMModel.Sigma(i,i,k));
                    pd = makedist('Normal');
                    pd.mu=my_mu;
                    pd.sigma=my_sig;
                    xmin=min(X_raw(:,i));
                    xmax=max(X_raw(:,i));
                    x = xmin:0.1:xmax;
                    pdf_normal = pdf(pd,x);
                    plot(x,pdf_normal*GMModel.ComponentProportion(k),'r','LineWidth',2);
                    sum_pdf=sum_pdf+pdf_normal*GMModel.ComponentProportion(k);
                    figure(num_fig)
                    subplot(sy,sy,(i-1)*sy+j)
                    
                end
                plot(x,sum_pdf,'k','LineWidth',2)
                histogram(X_raw(:,i),50,'Normalization','pdf')
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
                axis([min(X_raw(:,i)) max(X_raw(:,i)) 0 1.2])
                
                figure(num_fig)
                subplot(sy,sy,(i-1)*sy+j)
                
            else
                

                hold on
                [Z,~,~] = hist2(X_raw(:,i),X_raw(:,j),50,50);
                min1=min(X_raw(:,j));
                min2=min(X_raw(:,i));
                max1=max(X_raw(:,j));
                max2=max(X_raw(:,i));
                
                imagesc(Z)
                caxis([0 max(Z(:))-0.7*max(Z(:))])
                
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
                for k=1:nb_clust
                    my_mu1=GMModel.mu(k,j);
                    my_sig1=sqrt(GMModel.Sigma(j,j,k));
                    my_mu2=GMModel.mu(k,i);
                    my_sig2=sqrt(GMModel.Sigma(i,i,k));
                    my_cov=GMModel.Sigma(j,i,k);
                    Sig=[[my_sig1^2, my_cov];[my_cov, my_sig2^2]];
                    Mu=[my_mu1,my_mu2];
                    
                    x1 = min1:(max1-min1)/50:max1; x2 = min2:(max2-min2)/50:max2;
                    [X1,X2] = meshgrid(x1,x2);
                    F = mvnpdf([X1(:),X2(:)],Mu,Sig);
                    F = reshape(F,length(x2),length(x1));
                    
                    contour(F,[0.2 0.2001],'r')
                    
                end
                
            end
        end
    end
end

subplot(10,10,1)
xlabel('II.1')
subplot(10,10,12)
xlabel('II.2')
subplot(10,10,23)
xlabel('II.3')
subplot(10,10,34)
xlabel('SI.1')
subplot(10,10,45)
xlabel('SI.2')
subplot(10,10,56)
xlabel('SI.3')
subplot(10,10,67)
xlabel('SI.4')
subplot(10,10,78)
xlabel('TI.1')
subplot(10,10,89)
xlabel('TI.2')
subplot(10,10,100)
xlabel('TI.3')
end