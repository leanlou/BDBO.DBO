% 定义bdbo函数
function [Multi_Convergence_curve, Multi_time, bestfit] = bdbo(Z,Popsize, dim, ub, lb, fobj, Max_iteration, Multi_Num,func_num)
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
        %x = arrayfun(@(x) vpa(x, 20), x);
        fit = fobj(fhd,x',func_num);
        pFit = fit;
        pX = x;
        XX = pX;
        fMin = min(fit);
        bestI = find(fit == fMin, 1);
        bestX = x(bestI, :);
        unch_count = 0;
        Lambda = 1.5; 
        means = zeros(6, 1);
        covariances = ones(6, 6).*0.1;
        


        % 迭代过程
        for t = 1:Max_iteration
            % 最差的位置
            worse = x(find(fit == max(fit), 1), :);
            
            %滚球蜣螂更新
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
            R = PAO(t,Max_iteration); % S型曲线递减
            Xnew1 = max(min(bestX_ruc * (1 - R), ub), lb);
            Xnew2 = max(min(bestX_ruc * (1 + R), ub), lb);

            % 产卵区域
            Xnew11 = max(min(bestX * (1 - R), ub), lb);
            Xnew22 = max(min(bestX * (1 + R), ub), lb);

            % 育雏蜣螂更新
            

            if ii>Multi_Num*0.8 && t > 0.1*Max_iteration
                b1 = normrnd(means(1), sqrt(covariances(1,1)), [1, dim]);
                b2 = normrnd(means(2), sqrt(covariances(2,2)), [1, dim]);
            else
                b1=rand(1, dim);
                b2=rand(1, dim);
            end
            for i = rolNum+1:rolNum+babNum
                x(i, :) = bestX_ruc + b1 .* (pX(i, :) - Xnew1) + b2 .* (pX(i, :) - Xnew2);
                %x(i, :) = bestX + b1.* (pX(i, :) - Xnew1) + b2 * (pX(i, :) - Xnew2);
                x(i, :) = max(min(x(i, :), Xnew2), Xnew1);
                fit(i) = fobj(fhd,x(i, :)',func_num);
            end

            % 觅食蜣螂更新
            
            if ii>Multi_Num*0.8 && t > 0.1*Max_iteration
                C1 = normrnd(0, 1, [1, dim]);
                C2 = normrnd(means(3), sqrt(covariances(3,3)), [1, dim]);
            else
                C1=normrnd(0, 1, [1, dim]);
                C2=rand(1, dim);
            end
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
                    unch_count = 0;
                else
                    unch_count = unch_count + 1;%统计最优解未被更新次数
                end
            end
            
            

            % Update parameters
            %Levy飞行远距离更新
            if unch_count >=10 
                for i = 1:rolNum
                    % Scale factor
                    theta = 0.01 * (x(i,:) - bestX);
                    % Generate Levy step size
                    step_size = levy(Lambda, dim);
                    % Update position
                    x(i,:) = x(i,:) + theta .* step_size;
                    x(i, :) = max(min(x(i, :), ub), lb);
                end
                unch_count = 0;
            end

            Convergence_curve = [Convergence_curve; fMin];
            
        end

        % Return the results
        Multi_Convergence_curve = [Multi_Convergence_curve; Convergence_curve];
        Multi_time = [Multi_time; toc(start_time)];
        bestfit=[bestfit;fMin];
    end
end

%定义目标函数
function fit=fobj(fhd,x,func_num)
    %y = sum(x.^2);
    fit= feval(fhd, x,func_num);
end

% 定义抛物线型R函数
function pao = PAO(t,Max_iteration)
    v = sqrt(5.*Max_iteration);
    t0 = t/v;
    pao = 1- (0.5 .*10 .*(t0).^2)/Max_iteration;
end

% 定义Levy飞行函数
function step = levy(Lambda, dim)
    sigma = (gamma(1 + Lambda) * sin(pi * Lambda / 2) / (gamma((1 + Lambda) / 2) * Lambda * 2 ^ ((Lambda - 1) / 2))) ^ (1 / Lambda);
    u = randn(1, dim) * sigma;
    v = randn(1, dim);
    step = u ./ abs(v) .^ (1 / Lambda);
end
