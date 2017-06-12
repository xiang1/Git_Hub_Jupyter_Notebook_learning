% document： E:\嗑研专用文件夹\毫米波损耗研究\journal article\期刊相关文档\20170511备选的仿真.docx

%% 1.1 近似的效果仿真
%% 1.1.1 大尺度的效果 vs 阴影标准差
    %1.1.1.A  主要探索两个方面，其一是随着阴影标准差的变化，LoS在小区前半部分的近似效果，需要和LS方法进行对比
            clear; clc
            global glo_PARA
            glo_PARA.R = 400; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;
            i = 0;
            fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-10, 'TolX', 1e-10);
            loss = 60:0.1:180;
            sigma_2V = [5, 5.8, 8.7];
            FMatrix = zeros(length(loss), length(sigma_2V));
            FstFMatrix = zeros(length(loss), length(sigma_2V));
            ScdFMatrix = zeros(length(loss), length(sigma_2V));
            LSFMatrix = zeros(length(loss), length(sigma_2V));
            for i=1:1:length(sigma_2V)
                glo_PARA.sigma_2 = sigma_2V(i);
                [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();    %f21_2_l = f21_2_l.*sum(f21_l_pri)/sum(f21_2_l);
%                 [f22_l_pri,f22_2_l] = FunctionGetF22();    f22_2_l = f22_2_l.*sum(f22_l_pri)/sum(f22_2_l);
%                 [f31_l_pri,f31_2_l] = FunctionGetF31();    f31_2_l = f31_2_l.*sum(f31_l_pri)/sum(f31_2_l);
%                 [f32_l_pri,f32_2_l] = FunctionGetF32();    f32_2_l = f32_2_l.*sum(f32_l_pri)/sum(f32_2_l);

                f_pri = f21_l_pri./sum(f21_l_pri);%+f22_l_pri+f31_l_pri+f32_l_pri;
                f_1_l= f21_l./sum(f21_l);
                f_2_l = f21_2_l./sum(f21_2_l);%  +f22_2_l  +f31_2_l   +f32_2_l;
                FMatrix(:, i) = f_pri;
                FstFMatrix(:,i) = f_1_l';
                ScdFMatrix(:,i)= f_2_l';
                % LS fitting
                f_pri = f_pri./sum(f_pri);
                [fitmodel, gof] = fit(loss', f_pri, 'gauss1', fo);
                fG = fitmodel(loss');
                LSFMatrix(:,i) = fG';
            end;
            plot(loss, FMatrix(:,1), 'r',loss, ScdFMatrix(:,1), 'r--', loss, LSFMatrix(:,1), 'r:',...
                loss, FMatrix(:,2),'k',loss, ScdFMatrix(:,2),'k--',loss, LSFMatrix(:,2),'k:',...
                loss, FMatrix(:, 3),'b',loss, ScdFMatrix(:, 3), 'b--',loss, LSFMatrix(:, 3), 'b:', 'linewidth', 2.0);
            grid on; xlabel('loss(dB)'); ylabel('PDF');
            legend('Monte Carlo Simmulation \sigma_{S_2} = 5','2nd Step Approximation \sigma_{S_2} = 5',...
                'LS Fitting Approximation \sigma_{S_2} = 5',...
                'Monte Carlo Simmulation \sigma_{S_2} = 5.8','2nd Step Approximation \sigma_{S_2} = 5.8',...
                'LS Fitting Approximation \sigma_{S_2} = 5.8',...
                'Monte Carlo Simmulation \sigma_{S_2} = 8.7', '2nd Step Approximation \sigma_{S_2} = 8.7',...
                'LS Fitting Approximation \sigma_{S_2} = 8.7');
            KLV =zeros(3, length(sigma_2V));
            FMatrix = FMatrix./sum(FMatrix); FstFMatrix = FstFMatrix./sum(FstFMatrix); 
            ScdFMatrix= ScdFMatrix./sum(ScdFMatrix); LSFMatrix = LSFMatrix;%./sum(LSFMatrix);
            KLV(1,:) = -sum(FstFMatrix.*log(FMatrix./(FstFMatrix+eps)));
            KLV(2,:) = -sum(ScdFMatrix.*log(FMatrix./ScdFMatrix));
            KLV(3,:) = -sum(LSFMatrix.*log(FMatrix./LSFMatrix));
    %% 1.1.1 大尺度的效果 vs 路损因子
        %1.1.1.B  主要探索两个方面，其二是随着路损因子的变化，LoS在小区前半部分的近似效果，需要和LS方法进行对比
            clear; clc
            global glo_PARA
            glo_PARA.R = 400; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;
            i = 0;
            fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-10, 'TolX', 1e-10);
            loss = 60:0.1:180;
            beta_2V = [ 2.0, 2.5,3,3.5];
            FMatrix = zeros(length(loss), length(beta_2V));
            FstFMatrix = zeros(length(loss), length(beta_2V));
            ScdFMatrix = zeros(length(loss), length(beta_2V));
            LSFMatrix = zeros(length(loss), length(beta_2V));
            for i=1:1:length(beta_2V)
                glo_PARA.beta_2 = beta_2V(i);
                [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();    %f21_2_l = f21_2_l.*sum(f21_l_pri)/sum(f21_2_l);
%                 [f22_l_pri,f22_2_l] = FunctionGetF22();    f22_2_l = f22_2_l.*sum(f22_l_pri)/sum(f22_2_l);
%                 [f31_l_pri,f31_2_l] = FunctionGetF31();    f31_2_l = f31_2_l.*sum(f31_l_pri)/sum(f31_2_l);
%                 [f32_l_pri,f32_2_l] = FunctionGetF32();    f32_2_l = f32_2_l.*sum(f32_l_pri)/sum(f32_2_l);

                f_pri = f21_l_pri./sum(f21_l_pri);%+f22_l_pri+f31_l_pri+f32_l_pri;
                f_1_l= f21_l./sum(f21_l);
                f_2_l = f21_2_l./sum(f21_l);%  +f22_2_l  +f31_2_l   +f32_2_l;
                FMatrix(:, i) = f_pri;
                FstFMatrix(:,i) = f_1_l';
                ScdFMatrix(:,i)= f_2_l';
                % LS fitting
                f_pri = f_pri./sum(f_pri);
                [fitmodel, gof] = fit(loss', f_pri, 'gauss1', fo);
                fG = fitmodel(loss');
                LSFMatrix(:,i) = fG';
            end;
            plot(loss, FMatrix(:,1), 'r',loss, ScdFMatrix(:,1), 'r--', loss, LSFMatrix(:,1), 'r:',...
                loss, FMatrix(:,2),'k',loss, ScdFMatrix(:,2),'k--',loss, LSFMatrix(:,2),'k:',...
                loss, FMatrix(:, 3),'b',loss, ScdFMatrix(:, 3), 'b--',loss, LSFMatrix(:, 3), 'b:',... 
                loss, FMatrix(:, 4),'m',loss, ScdFMatrix(:, 4), 'm--',loss, LSFMatrix(:, 4), 'm:' ,'linewidth', 2.0);
            grid on; xlabel('loss(dB)'); ylabel('PDF');
            legend('Monte Carlo Simmulation \beta_2 = 2','2nd Step Approximation \beta_2 = 2',...
                'LS Fitting Approximation \beta_2 = 2',...
                'Monte Carlo Simmulation \beta_2 = 2.5','2nd Step Approximation \beta_2 = 2.5',...
                'LS Fitting Approximation \beta_2 = 2.5',...
                'Monte Carlo Simmulation \beta_2 = 3', '2nd Step Approximation \beta_2 = 3',...
                'LS Fitting Approximation \beta_2 = 3',...
                'Monte Carlo Simmulation \beta_2 = 3.5', '2nd Step Approximation \beta_2 = 3.5',...
                'LS Fitting Approximation \beta_2 = 3.5');
            KLV =zeros(3, length(beta_2V));
            FMatrix = FMatrix./sum(FMatrix); FstFMatrix = FstFMatrix./sum(FstFMatrix); 
            ScdFMatrix= ScdFMatrix./sum(ScdFMatrix); LSFMatrix = LSFMatrix./sum(LSFMatrix);
            KLV(1,:) = -sum(FstFMatrix.*log(FMatrix./(FstFMatrix+eps)));
            KLV(2,:) = -sum(ScdFMatrix.*log(FMatrix./ScdFMatrix));
            KLV(3,:) = -sum(LSFMatrix.*log(FMatrix./LSFMatrix));
            KLV = KLV';
	%% 1.1.2 大尺度+小尺度损耗的pdf
        % 考虑用Nagakami-m的统一信道模型，
        % 改变m值即可认为是不同的信道m=1 瑞利信道，m = 1/(1-(K/(1+K))^2 ) 认为是Rice信道
        % 根据相关文献，los下莱斯分布的k值范围未7dB~17dB
            clear; clc
            fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-10, 'TolX', 1e-10);
            global glo_PARA
            glo_PARA.R = 400; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;
            KV = [7, 12, 17];
            global loss
            loss = 60:0.1:180 ;
            FMatrixSmFd = zeros(length(loss), length(KV));
            ScdFMatrixSmFd = zeros(length(loss), length(KV));
            LSMatrixSFd = zeros(length(loss), length(KV));
            beta_2V = [2,2.5,3];
            for i = 1:1:length(KV)
                K = KV(i);
                glo_PARA.beta_2 = beta_2V(i); 
                [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();  
                pr = sum(f21_l_pri);
                f21_l_pri = f21_l_pri./sum(f21_l_pri);
                f21_2_l = f21_2_l./sum( f21_2_l );
                [fLoss, fApprox] = FunctionLossWithSmFd(f21_l_pri,f21_2_l, K );
                fLoss = fLoss ./sum(fLoss);
                fApprox = fApprox./sum(fApprox);
                FMatrixSmFd(:,i) = fLoss';
                ScdFMatrixSmFd(:,i) = fApprox';
                [fitmodel, gof] = fit(loss', fLoss, 'gauss1', fo);
                fG = fitmodel(loss');
                LSMatrixSFd(:,i) = fG';
            end
            
            plot( loss,FMatrixSmFd(:,1),'r', loss,  ScdFMatrixSmFd(:,1),'r--',loss, LSMatrixSFd(:,1),'r:',...
                loss,FMatrixSmFd(:,2),'b', loss,  ScdFMatrixSmFd(:,2),'b--',loss, LSMatrixSFd(:,2),'b:',...
                loss,FMatrixSmFd(:,3),'k', loss,  ScdFMatrixSmFd(:,3),'k--',loss, LSMatrixSFd(:,3),'k:','linewidth',1.5);
            grid on;  xlabel('loss(dB)'); ylabel('PDF');
            legend('Disttribution of L_{21} with k =7, \beta_2 = 2',' Approximation with k =7,\beta_2 = 2','LS Fitting with k =7,\beta_2 = 2',...
                        'Disttribution of L_{21} with k =12, \beta_2 = 2.5',' Approximation with k =12, \beta_2 = 2.5','LS Fitting with k =12, \beta_2 = 2.5',...
                        'Disttribution of L_{21} with k =17, \beta_2 = 3',' Approximation with k =17, \beta_2 = 3','LS Fitting with k =12, \beta_2 = 2.5')
            KLV(1,:) = -sum(ScdFMatrixSmFd.*log(FMatrixSmFd./(ScdFMatrixSmFd+eps)));
            KLV(2,:) = -sum(LSMatrixSFd.*log(FMatrixSmFd./LSMatrixSFd));        
%% 1.2 统计量的仿真 均值 & 方差 vs 小区半径
       % 仿真不同状态下的损耗均值和方差
            clear; clc
            global glo_PARA
            glo_PARA.R = 400; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;
            RV = 50:10:400;
            loss = 60:0.1:180;
            meanMatrix = zeros( length(RV),6);
            varMatrix = zeros(  length(RV),6);
            for i = 1: 1: length(RV)
                 glo_PARA.R = RV(i);
                 K_los = 12; K_nlos  = 15;
                 [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();            [f21_l_pri, f21_2_l] = FunctionLossWithSmFd(f21_l_pri,f21_2_l, K_los );
                [f22_l_pri,f22_2_l] = FunctionGetF22();                       [f22_l_pri, f22_2_l] = FunctionLossWithSmFd(f22_l_pri,f22_2_l, K_los );
                [f31_l_pri,f31_2_l] = FunctionGetF31();                       [f31_l_pri, f31_2_l] = FunctionLossWithSmFd(f31_l_pri,f31_2_l, K_nlos );
                [f32_l_pri,f32_2_l] = FunctionGetF32();                       [f32_l_pri, f32_2_l] = FunctionLossWithSmFd(f32_l_pri,f32_2_l, K_nlos );
                
                f_pri = f21_l_pri+f22_l_pri+f31_l_pri+f32_l_pri;
                f_2_l = f21_2_l+f22_2_l  +f31_2_l   +f32_2_l;
                meanMatrix(i,1 ) = sum(loss.* f_pri' )*0.1;%./sum(f_pri));
                meanMatrix(i,2 ) = sum(loss.* f_2_l )*0.1;%./sum(f_2_l));
                meanMatrix(i,3 ) = sum(loss.* (f21_l_pri'+f22_l_pri')./sum((f21_l_pri'+f22_l_pri')));
                meanMatrix(i,4 ) = sum(loss.* (f21_2_l+f22_2_l)./sum( (f21_2_l+f22_2_l)));
                meanMatrix(i,5 ) = sum(loss.* (f31_l_pri'+f32_l_pri')./sum((f31_l_pri'+f32_l_pri')));
                meanMatrix(i,6 ) = sum(loss.* (f31_2_l+f32_2_l)./sum(f31_2_l+f32_2_l));

%                 varVec = (repmat(loss, 6,1)- repmat(meanMatrix(i,:)',1, length(loss))).^2;
%                 varMatrix(i,1 ) = sum(varVec(1,:).* f_pri' )*0.1;%./sum(f_pri));
%                 varMatrix(i,2 ) = sum(varVec(2,:).* f_2_l )*0.1;%./sum(f_2_l));
%                 varMatrix(i,3 ) = sum(varVec(3,:).* (f21_l_pri'+f22_l_pri')./sum((f21_l_pri'+f22_l_pri')));
%                 varMatrix(i,4 ) = sum(varVec(4,:).* (f21_2_l+f22_2_l)./sum( (f21_2_l+f22_2_l)));
%                 varMatrix(i,5 ) = sum(varVec(5,:).* (f31_l_pri'+f32_l_pri')./sum(f31_l_pri'+f32_l_pri'));
%                 varMatrix(i,6 ) = sum(varVec(6,:).* (f31_2_l+f32_2_l)./sum(f31_2_l+f32_2_l));
                 varVec = (repmat(loss, 6,1)).^2;
                varMatrix(i,1 ) = sum(varVec(1,:).* f_pri' )*0.1;%./sum(f_pri));
                varMatrix(i,2 ) = sum(varVec(2,:).* f_2_l  )*0.1;%./sum(f_2_l));
                varMatrix(i,3 ) = sum(varVec(1,:).* (f21_l_pri'+f22_l_pri')  ./sum((f21_l_pri'+f22_l_pri')));
                varMatrix(i,4 ) = sum(varVec(2,:).* (f21_2_l+f22_2_l)   ./sum( (f21_2_l+f22_2_l)));
                varMatrix(i,5 ) = sum(varVec(1,:).* (f31_l_pri'+f32_l_pri')     ./sum(f31_l_pri'+f32_l_pri'));
                varMatrix(i,6 ) = sum(varVec(2,:).* (f31_2_l+f32_2_l)   ./sum(f31_2_l+f32_2_l));
%       
                
            end
            plot(RV,meanMatrix(:,1),'ok',...
                    RV,meanMatrix(:,2),'k',...
                    RV,meanMatrix(:,3),'rx',...
                    RV,meanMatrix(:,4),'r--',...
                    RV,meanMatrix(:,5),'bp',...
                    RV,meanMatrix(:,6),'b:',...        
                     'linewidth',2.0);grid minor;
            xlabel('Radius of the cell R(m)');ylabel('First moment(mean) of loss(dB)');
            legend('System','Approximation of the system',...
                'LoS contribution', 'Approximation of the LoS state',...
                'NLoS contribution', 'Approximation of the NLoS state');
            figure
            %varMatrix = sqrt(varMatrix);
            plot(RV,varMatrix(:,1),'ok',...
                    RV,varMatrix(:,2),'k',...
                    RV,varMatrix(:,3),'rx',...
                    RV,varMatrix(:,4),'r--',...
                    RV,varMatrix(:,5),'bp',...
                    RV,varMatrix(:,6),'b:',...        
                     'linewidth',2.0);grid minor;
            xlabel('Radius of the cell R(m)');ylabel('Second moment of loss');

            legend('System','Approximation of the system',...
                'LoS contribution', 'Approximation of the LoS state',...
                'NLoS contribution', 'Approximation of the NLoS state');
%% 1.3 系统性能的仿真
    %% 1.3.1 a 覆盖率 vs 阈值 for different 小区半径
            %覆盖率在不同半径条件下随着阈值的变化情况
            clc;clear;
            global glo_PARA
            glo_PARA.R = 200; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;

            %subplot(2,2,1);
            i = 0;
            rate_pri = zeros(1,11);
            rate_2=zeros(1,11);
            RV = [50.0,100,200];
            original=zeros(length(1:50:length(60:0.1:180)),length(200:40:350));
            approx = zeros(length(60:0.1:180),length(200:40:400));
             legendStr=[];
             legendStr2=[];
            for i = 1:1:length(RV)
                glo_PARA.R  = RV(i);
                 K_los = 12; K_nlos  = 15;
                loss = 60:0.1:180;
                plotloss =1:50:length(loss);
                 [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();            [f21_l_pri, f21_2_l] = FunctionLossWithSmFd(f21_l_pri,f21_2_l, K_los );
                [f22_l_pri,f22_2_l] = FunctionGetF22();                       [f22_l_pri, f22_2_l] = FunctionLossWithSmFd(f22_l_pri,f22_2_l, K_los );
                [f31_l_pri,f31_2_l] = FunctionGetF31();                       [f31_l_pri, f31_2_l] = FunctionLossWithSmFd(f31_l_pri,f31_2_l, K_nlos );
                [f32_l_pri,f32_2_l] = FunctionGetF32();                       [f32_l_pri, f32_2_l] = FunctionLossWithSmFd(f32_l_pri,f32_2_l, K_nlos );

                f_pri = f21_l_pri+f22_l_pri+f31_l_pri+f32_l_pri;
                f_1_l = f21_2_l+f22_2_l  +f31_2_l   +f32_2_l;
                plotf_pri = cumsum(f_pri);%+f_pri(end);
                original(:,i) = plotf_pri(plotloss)*0.1;
                approx(:,i) = cumsum(f_1_l)*0.1;
            end;
            originalPlot = original;
            approxPlot = approx;
            SNR = -loss+30+20+70.9897;    %30是发射天线 的功率 20 天线增益 160.89噪声功率 
            H1 =plot(SNR(plotloss),originalPlot(:,1),'ko',...
                          SNR(plotloss),originalPlot(:,2),'r pentagram',...
                          SNR(plotloss),originalPlot(:,3),'bx','linewidth',2.0);hold on;

            xlabel('SNR Threshold(dB)');ylabel('P_{cov}'); 
            set(gca,'FontName','Helvetica')

            H2  = plot(SNR,approxPlot(:,1),'k',...
                            SNR,approxPlot(:,2),'r--',...
                            SNR,approxPlot(:,3),'b:','linewidth',2.0);
            hold on;            grid on;
            legend('R=50m,Monte Carlo Simulation',...
                        'R=100m,Monte Carlo Simulation',...
                        'R=200m,Monte Carlo Simulation',...
                        'R=50m,Approximation',...
                        'R=100m,Approximation',...
                        'R=200m,Approximation');
    %% 1.3.1 b  覆盖率 vs 小区半径 for 特定的阈值
            %给定覆盖的阈值，不同状态下对覆盖率的贡献变化曲线
            clc;clear;
            global glo_PARA
            glo_PARA.R = 200; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;

            figure;
            i = 0;
            RV = 50:15:400;
            outageRate = zeros(length(RV),6);
            coverageWinnerV= zeros(length(RV),1);
            threshold =-20.00103;
            for i = 1:1:length(RV)
                glo_PARA.R = RV(i);
                 K_los = 12; K_nlos  = 15;
                loss = 60:0.1:180;
                plotloss =1:50:length(loss);
                [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();            [f21_l_pri, f21_2_l] = FunctionLossWithSmFd(f21_l_pri,f21_2_l, K_los );
                [f22_l_pri,f22_2_l] = FunctionGetF22();                       [f22_l_pri, f22_2_l] = FunctionLossWithSmFd(f22_l_pri,f22_2_l, K_los );
                [f31_l_pri,f31_2_l] = FunctionGetF31();                       [f31_l_pri, f31_2_l] = FunctionLossWithSmFd(f31_l_pri,f31_2_l, K_nlos );
                [f32_l_pri,f32_2_l] = FunctionGetF32();                       [f32_l_pri, f32_2_l] = FunctionLossWithSmFd(f32_l_pri,f32_2_l, K_nlos );

                f_pri = f21_l_pri+f22_l_pri+f31_l_pri+f32_l_pri;
                f_1_l = f21_2_l+f22_2_l  +f31_2_l   +f32_2_l;
                loss = 60:0.1:180;
                SNR = -loss+30+20+70.9897;
                outageRate(i,1) = 0.1*sum(f_pri(1:max(find(SNR>threshold))));
                outageRate(i,2) = 0.1*sum(f_1_l(1:max(find(SNR>threshold))));
                f_tmp =  f21_l_pri+f22_l_pri;
                outageRate(i,3) = 0.1*sum(f_tmp(1:max(find(SNR>threshold))));
                f_tmp = f21_2_l+f22_2_l;
                outageRate(i,4) = 0.1*sum(f_tmp(1:max(find(SNR>threshold))));
                f_tmp =  f31_l_pri+f32_l_pri;
                outageRate(i,5) = 0.1*sum(f_tmp(1:max(find(SNR>threshold))));
                    f_tmp =f31_2_l   +f32_2_l;
                outageRate(i,6) = 0.1*sum(f_tmp(1:max(find(SNR>threshold))));
            end;
            H3 = plot(RV, outageRate(:,1),'ko',...
                    RV, outageRate(:,2),'k',...
                    RV, outageRate(:,3),'rx',...
                    RV, outageRate(:,4),'r--',...
                    RV, outageRate(:,5),'bpentagram',...
                    RV, outageRate(:,6),'b:',...
                    'linewidth',2.0);
           grid on;
            xlabel('Radius of the cell R(m)');ylabel('Coverage');
            legend('Coverage of the system', 'Approximation for the coverage of the system',...
                            'Coverage of the LoS contribution','Approximation for the LoS contribution',...
                            'Coverage of the NLoS contribution', 'Approximation for the NLoS contribution');
            figure;
            plot(RV, outageRate(:,3)./outageRate(:,1)*100,'rx',...
                RV, outageRate(:,4)./outageRate(:,2)*100,'r--',...
                RV, outageRate(:,5)./outageRate(:,1)*100,'bp',...
                RV,outageRate(:,6)./outageRate(:,2)*100,'b:','linewidth', 2.0);
            xlabel('radius of the cell (m)'); ylabel('contributions(%)');
            grid on;
            legend('The fraction of LoS coverage contribution',...
                'Approximation for the fraction LoS contribution',...
                'The fraction of NLoS coverage contribution',...
                'Approximation for the fraction NLoS contribution')
%% 1.3.2 a 香农公式的近似
       % log_2(1+x) -> a+b*log10(x) +c*( log10(x) )^2
           clear;clc;
           x= 0.0001:0.1:100;
           fx = log2(1+x);
           fo = fitoptions('Method', 'NonlinearLeastSquares',...
                            'MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-10, 'TolX', 1e-10);
            [fitmodel, gof] = fit(x', fx', 'a+b*log10(x) +c*( log10(x) )^2', fo);
            fG = fitmodel(x);
            %fG = 1.153 + 1.843*log10(x) +0.4709*(log10(x)).^2;
            plot(x, fx,'k',x, fG,'r--','linewidth',2.0);
            grid on; xlabel('x'); ylabel('y')
            legend(' y = log_2(1+x)',' y = a_0+a_1 log_{10}(x) +a_2 (log_{10}(x) )^2');
%% 1.3.2 b 微波毫米波 rate cdf对比 vs 小区半径
        %%
        clc;clear;
        global glo_PARA
        glo_PARA.R = 200; 
        glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
        glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
        glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;
        
        figure;        i = 0;
        rate_pri = zeros(1,11);
        rate_2=zeros(1,11);
        RV = [50,100,200];%,300,400];
        original=zeros(length(RV),length(60:0.1:180));
        approx = zeros(length(RV),length(60:0.1:180));
        original_Mu_Wave = zeros(length(RV), length(20:0.1:180));
         legendStr=[];
         legendStr2=[];
         rate_Results = zeros(length(RV)*2, 3);
        for i = 1:1:length(RV)
            glo_PARA.R  = RV(i);
            % R = 400;
            loss = 60:0.1:180;
            plotloss =1:50:length(loss);
            K_los = 12; K_nlos  = 15;
            [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();             [f21_l_pri, f21_2_l] = FunctionLossWithSmFd(f21_l_pri,f21_2_l, K_los );
            [f22_l_pri,f22_2_l] = FunctionGetF22();                       [f22_l_pri, f22_2_l] = FunctionLossWithSmFd(f22_l_pri,f22_2_l, K_los );
            [f31_l_pri,f31_2_l] = FunctionGetF31();                       [f31_l_pri, f31_2_l] = FunctionLossWithSmFd(f31_l_pri,f31_2_l, K_nlos );
            [f32_l_pri,f32_2_l] = FunctionGetF32();                       [f32_l_pri, f32_2_l] = FunctionLossWithSmFd(f32_l_pri,f32_2_l, K_nlos );

             f_pri = f21_l_pri+f22_l_pri+f31_l_pri+f32_l_pri;
             f_2_l = f21_2_l+f22_2_l  +f31_2_l   +f32_2_l;
            if(glo_PARA.R<156)
                f_pri = f_pri./sum(f_pri)*10;
                f_2_l = f_2_l./sum(f_2_l).*sum(f_pri);
            end
             Bw = 1000;
             SNR = -loss+30+20+174-12-10*log10(Bw*10^6);    %30是发射天线 的功率, 20 天线增益10, 接受天线的功率 160.89噪声功率
            rate =Bw.*(log2(1+10.^((SNR)./10)));
            % 微波
             loss2 = 20:0.1:180;    Bw2 = 20;
             SNR2 = -loss2+30+174-12-10*log10(Bw2*10^6);  % 微波 
            rate2 = Bw2.*(log2(1+10.^((SNR2)./10)));
            f_Mu_Wave = WinnerIIC2(glo_PARA.R );
            original_Mu_Wave(i,:) = 1-cumsum(f_Mu_Wave')*0.1;

            original(i,:) = 1-cumsum(f_pri')*0.1;
            approx(i,:) = 1-cumsum(f_2_l)*0.1;
            % 毫米波和微波的比较结果
            rate_Results(i,1) = sum(rate.*f_pri').*0.1/Bw;
            rate_Results(i,2) = sum(rate.*f_pri').*0.1;       
            rate_Results(i,3) = rate(find(original(i,:)>0.05,1,'last'));  %边缘用户

            rate_Results(i+length(RV),1) = sum(rate2.*f_Mu_Wave')*0.1/Bw2;
            rate_Results(i+length(RV),2) = sum(rate2.*f_Mu_Wave')*0.1;       
            rate_Results(i+length(RV),3) = rate2(find(original_Mu_Wave(i,:)>0.05,1,'last'));  %边缘用户
        end;
        ratePlot= 1:50:length(rate);

        H1 =semilogx(rate(ratePlot),original(1,ratePlot),'ko',...
                      rate(ratePlot),original(2,ratePlot),'rp',...                      %rate(ratePlot),original(3,ratePlot),'bx',...
                      rate,approx(1,:),'k',...
                      rate,approx(2,:),'r--',...                      %rate,approx(3,:),'b:',...
                      rate2,original_Mu_Wave(1,:),...
                      rate2, original_Mu_Wave(2,:),...                      %rate2, original_Mu_Wave(3,:),...
                      'linewidth',2.0);hold on;

        xlim([0, 50000]);
        ylim([-.01,1])
        grid minor;
        xlabel('rate( Mbps)');ylabel('CDF');
        %title('CDF of rate for different Radius R');         
        legend('R = 50m, Monte Carlo Simulation','R = 100m,Monte Carlo Simulation','R = 200m, Monte Carlo Simulation',...
                    'R = 50m, Approximation','R = 100m, Approximation','R = 200m, Approximation');


%% 1.3.2.c 速率 vs 小区半径
            clear; clc;
            % log_2(1+x) \approx  1.137+1.868*log10(x) +0.4621*( log10(x) )^2
            global glo_PARA
            glo_PARA.R = 200; 
            glo_PARA.alpha_los = 1/67.1;    glo_PARA.alpha_out = 1/30;  glo_PARA.beta_out = 5.2;
            glo_PARA.alpha_2 = 61.4;          glo_PARA.beta_2 = 2;            glo_PARA.sigma_2 = 5.8;
            glo_PARA.alpha_3 = 72;            glo_PARA.beta_3 = 2.92;        glo_PARA.sigma_3 = 8.7;

            P =30;% dbm发射功率
            Ndb = 10;% db
            RV = 50:15:400;
            rate_pri = zeros(size(RV));
            rate_2=zeros(size(RV));
            rate_los = zeros(size(RV));
            rateA_los = zeros(size(RV));
            rate_nlos = zeros(size(RV));
            rateA_nlos = zeros(size(RV));
            for i = 1:1:length(RV)
                glo_PARA.R = RV(i);
            % R = 400;
                loss = 60:0.1:180;
                K_los = 12; K_nlos  = 15;
                [f21_l_pri, f21_l, f21_2_l] = FunctionGetF21();            [f21_l_pri, f21_2_l] = FunctionLossWithSmFd(f21_l_pri,f21_2_l, K_los );
                [f22_l_pri,f22_2_l] = FunctionGetF22();                       [f22_l_pri, f22_2_l] = FunctionLossWithSmFd(f22_l_pri,f22_2_l, K_los );
                [f31_l_pri,f31_2_l] = FunctionGetF31();                       [f31_l_pri, f31_2_l] = FunctionLossWithSmFd(f31_l_pri,f31_2_l, K_nlos );
                [f32_l_pri,f32_2_l] = FunctionGetF32();                       [f32_l_pri, f32_2_l] = FunctionLossWithSmFd(f32_l_pri,f32_2_l, K_nlos );


                f_pri = f21_l_pri+f22_l_pri+f31_l_pri+f32_l_pri;
                f_2_l = f21_2_l+f22_2_l  +f31_2_l   +f32_2_l;
                
                SNR = -loss+30+20+70.9897;
                SNR_line = 10.^((SNR)./10);
                Rate_of_SNR = log2(1+SNR_line);
                Rate_of_SNR_App =Rate_of_SNR;% 1.137+1.868*log10(SNR_line) +0.4621*( log10(SNR_line) ).^2;
                
                rate_pri(i) = sum(Rate_of_SNR.*f_pri').*0.1;
                rate_2(i) = sum(Rate_of_SNR_App.*f_2_l).*0.1;

                rate_los(i) = sum(Rate_of_SNR.*(f21_l_pri'+f22_l_pri')).*0.1;
                rateA_los(i) = sum(Rate_of_SNR_App.*(f21_2_l+f22_2_l)).*0.1;

                rate_nlos(i) = sum(Rate_of_SNR.*(f31_l_pri'+f32_l_pri')).*0.1;
                rateA_nlos(i) = sum(Rate_of_SNR_App.*(f31_2_l+f32_2_l)).*0.1;

            end;
            %subplot(2,2,4)
            plot(RV(1:1:length(RV)),rate_pri,'ok',...
                    RV,rate_2,'k',...
                    RV, rate_los,'rx',...
                    RV,rateA_los,'r--',...
                    RV,rate_nlos,'bp',...
                    RV,rateA_nlos,'b:',...        
                     'linewidth',2.0);grid minor;
            xlabel('Radius of the cell R(m)');ylabel('Rate/Bw [bps/Hz]');
            % title('Average Rate vs Radius of cell');
            % set(gca,'FontName',)
            legend('Average rate of the system','Approximation for the average rate of the system',...
                'Average rate of LoS contribution', 'Approximation for the LoS contribution',...
                'Average rate of NLoS contribution', 'Approximation for the NLoS contribution');
           figure;
           plot(RV, rate_los./rate_pri*100,'rx',...
                RV, rateA_los./rate_2*100,'r--',...
                RV, rate_nlos./rate_pri*100,'bp',...
                RV,rateA_nlos./rate_2*100,'b:','linewidth', 2.0);
            xlabel('radius of the cell (m)'); ylabel('contributions(%)');
            grid on;
            legend('The fraction of LoS rate contribution',...
                'Approximation for the fraction LoS contribution',...
                'The fraction of NLoS rate contribution',...
                'Approximation for the fraction NLoS contribution')  
           