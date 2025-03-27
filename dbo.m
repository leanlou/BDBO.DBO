% Call the function
function [Multi_Convergence_curve, Multi_time, bestfit] = dbo(Z,Popsize, dim, ub, lb, fobj, Max_iteration, Multi_Num,func_num)
    Multi_Convergence_curve = []; 
    Multi_time = [];
    bestfit=[];

    for ii = 1:Multi_Num
        start_time = tic;
        % 定义四种蜣螂的数量
        rolNum = round(0.25 * Popsize); 
        eggNum = round(0.25 * Popsize);
        babNum = round(0.25 * Popsize);
        steNum = round(0.25 * Popsize);

        % 初始化种群
        fhd=str2func('cec17_func');
        seed = Z; % 设置你想要的种子值
        rng(seed, 'twister'); % 使用'Twister'算法设置随机数生成器的种子为seed 

        %x = rand(Popsize, dim) * (ub - lb) + lb;
        for i = 1 : Popsize
            x( i, : ) = lb + (ub - lb) .* rand( 1, dim );  
            %fit( i ) = fobj( x( i, : ) ) ;                       
        end
        rng('default')
        Convergence_curve = [];
        Convergence_curve=[Convergence_curve;min(fobj(fhd,x',func_num))];
        fit = fobj(fhd,x',func_num);
        pFit = fit;
        pX = x;
        XX = pX;
        fMin = min(fit);
        bestI = find(fit == fMin, 1);
        bestX = x(bestI, :);
        

        % 迭代过程
        for t = 1:Max_iteration
            % 最差的位置
            worse = x(find(fit == max(fit), 1), :);
            
            % 跳舞or滚球概率
            r2 = rand();
            for i = 1:rolNum
                if r2 < 0.9
                    % 滚球
                    r1 = rand();
                    % 自然系数
                    a = rand();
                    if a > 0.1
                        a = 1;
                    else
                        a = -1;
                    end
                    % 滚球蜣螂位置
                    x(i, :) = pX(i, :) + 0.3 * abs(pX(i, :) - worse) + a * 0.1 * XX(i, :);
                else
                    % 跳舞
                    aaa = randi([0, 180]);
                    if ismember(aaa, [0, 90, 180])
                        x(i, :) = pX(i, :);
                    end
                    theta = aaa * pi / 180;
                    x(i, :) = pX(i, :) + tan(theta) * abs(pX(i, :) - XX(i, :));
                end
                
                % 更新
                x(i, :) = max(min(x(i, :), ub), lb);
                fit(i) = fobj(fhd,x(i, :)',func_num);
            end
            

            %当前最优解
            fMin_ruc = min(fit);   % fMin_ruc是当前的最优适应度值
            bestII = find(fit == fMin_ruc, 1);
            bestX_ruc = x(bestII, :);% bestX_ruc是当前的最优解


           % 育雏球
            R = 1 - t/Max_iteration; % S型曲线递减
            Xnew1 = max(min(bestX_ruc * (1 - R), ub), lb);
            Xnew2 = max(min(bestX_ruc * (1 + R), ub), lb);

            % 产卵区域
            Xnew11 = max(min(bestX * (1 - R), ub), lb);
            Xnew22 = max(min(bestX * (1 + R), ub), lb);

            % 小蜣螂更新
            b1=rand(1, dim);
            b2=rand(1, dim);
            
            for i = rolNum+1:rolNum+babNum
                x(i, :) = bestX_ruc + b1 .* (pX(i, :) - Xnew1) + b2 .* (pX(i, :) - Xnew2);
                x(i, :) = max(min(x(i, :), Xnew2), Xnew1);
                fit(i) = fobj(fhd,x(i, :)',func_num);
            end

            % 育雏球更新
            
            
            C1=normrnd(0, 1, [1, dim]);
            C2=rand(1, dim);
            
            for i = rolNum+babNum+1:rolNum+babNum+eggNum
                x(i, :) = pX(i, :) + C1 .* (pX(i, :) - Xnew11) + C2 .* (pX(i, :) - Xnew22);
                x(i, :) = max(min(x(i, :), ub), lb);
                fit(i) = fobj(fhd,x(i, :)',func_num);
            end

            % 小偷蜣螂更新
            g=normrnd(0, 1, [1, dim]);
            
            for i = rolNum+babNum+eggNum+1:Popsize
                x(i, :) = bestX + g .* (abs(pX(i, :) - bestX_ruc) + abs(pX(i, :) - bestX)) / 2;
                x(i, :) = max(min(x(i, :), ub), lb);
                fit(i) = fobj(fhd,x(i, :)',func_num);
            end

            % 更新个体最优和全局最优
            for i = 1:Popsize
                if fit(i) < pFit(i)
                    pFit(i) = fit(i);
                    pX(i, :) = x(i, :);
                end
                if pFit(i) < fMin
                    fMin = pFit(i);
                    bestX = pX(i, :);
                end
            end
            Convergence_curve = [Convergence_curve;fMin];
        end

        % Return the results
        bestfit=[bestfit;fMin];
        Multi_Convergence_curve = [Multi_Convergence_curve; Convergence_curve];
        Multi_time = [Multi_time; toc(start_time)];
    end
end

% Call the function