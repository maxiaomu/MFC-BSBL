function Result = MFC_BSBL(Phi, y, blkStartLoc, LearnLambda, varargin)

% scaling...
scl = mean(std(y)); %std(y):每一列(默认)/行的标准偏差。得到一个向量再做均值
if (scl < 0.4) | (scl > 1)
    y = y/scl*0.4;
end

% get problem dimension
[M,N] = size(Phi);
L = size(y,2);

% Default Parameter Values for Any Cases
EPSILON       = 1e-6;       % solution accurancy tolerance  
MAX_ITERS     = 50;        % maximum iterations
PRINT         = 1;          % show progress information
LEARNTYPE     = 1;          % adaptively estimate the covariance matrix B

if LearnLambda == 0  
    lambda = 1e-12;   
    PRUNE_GAMMA = 1e-3;  MatrixReg = zeros(L);
elseif LearnLambda == 2
    lambda = 1e-3;    
    PRUNE_GAMMA = 1e-2;   MatrixReg = eye(L) * 2 ; % 1e-2变-3
elseif LearnLambda == 1
    lambda = 1e-3;    
    PRUNE_GAMMA = 1e-2;  MatrixReg = eye(L) * 4; %4
else
    error(['Unrecognized Value for LearnLambda']);
end

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'learntype'
                LEARNTYPE = varargin{i+1};
                if LEARNTYPE ~= 1 & LEARNTYPE ~= 0
                    error(['Unrecognized Value for LEARNTYPE']);
                end
            case 'prune_gamma'
                PRUNE_GAMMA = varargin{i+1}; 
            case 'lambda'
                lambda = varargin{i+1};    
            case 'epsilon'   
                EPSILON = varargin{i+1}; 
            case 'print'    
                PRINT = varargin{i+1}; 
            case 'max_iters'
                MAX_ITERS = varargin{i+1};  
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

if PRINT
    fprintf('\n====================================================\n');
    fprintf('           Running MFC-BSBL....... \n');
    fprintf('           Information about parameters...\n');
    fprintf('====================================================\n');
    fprintf('PRUNE_GAMMA  : %e\n',PRUNE_GAMMA);
    fprintf('lambda       : %e\n',lambda);
    fprintf('LearnLambda  : %d\n',LearnLambda);    
    fprintf('LearnType    : %d\n',LEARNTYPE);
    fprintf('EPSILON      : %e\n',EPSILON);
    fprintf('MAX_ITERS    : %d\n\n',MAX_ITERS);
end

count_num = 0; dmu = 0;
%% Initialization

y0 = y;
Phi0 = Phi;
blkStartLoc0 = blkStartLoc;
p = length(blkStartLoc);   % block number
for k = 1 : p-1
    blkLenList(k) = blkStartLoc(k+1)-blkStartLoc(k);
end
blkLenList(p) = N - blkStartLoc(end)+1;
maxLen = max(blkLenList);
if sum(blkLenList == maxLen) == p, 
    equalSize = 1;
else
    equalSize = 0;
end

for k = 1 : p
    Sigma0{k} = eye(blkLenList(k));
end

gamma = ones(p,1);
keep_list = [1:p]';
usedNum = length(keep_list);
mu_t = zeros(N,L);
count = 0;

B = eye(L);

%% Iteration
while (1)
    count = count + 1;

    %=========== Prune sparsity ==============
    if (min(gamma) < PRUNE_GAMMA)
        index = find(gamma > PRUNE_GAMMA);
        usedNum = length(index);
        keep_list = keep_list(index);
        if isempty(keep_list), 
            fprintf('\n====================================================================================\n');
            fprintf('x becomes zero vector. The solution may be incorrect. \n');
            fprintf('Current ''prune_gamma'' = %g, and Current ''EPSILON'' = %g.\n',PRUNE_GAMMA,EPSILON);
            fprintf('Try smaller values of ''prune_gamma'' and ''EPSILON'' or normalize ''y'' to unit norm.\n');
            fprintf('====================================================================================\n\n');
            break; 
        end;
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        
        % prune gamma and associated components in Sigma0 
        gamma = gamma(index);
        temp = Sigma0;
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        
        % construct new Phi
        temp = [];
        for k = 1 : usedNum
            temp = [temp, Phi0(:,blkStartLoc(k):blkStartLoc(k)+blkLenList(k)-1)];
        end
        Phi = temp;
    end

    % =================== Alternating-learning method ==============
    y = y0 * sqrtm(inv(B));
    
    %=================== Compute new weights =================
    PhiAPhi = zeros(M);
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        PhiAPhi = PhiAPhi + Phi(:, currentSeg)*Sigma0{i}*Phi(:, currentSeg)';
        currentLoc = currentSeg(end);
    end

    H = Phi' /(PhiAPhi + lambda * eye(M));
    Hy = H * y; 
    HPhi = H * Phi;
    
    mu_x = zeros(size(Phi,2),L);
    Sigma_x = [];
    
    A = []; invA = []; A0 = zeros(maxLen); r0 = zeros(1); r1 = zeros(1);
   
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        seg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        mu_x(seg,:) = Sigma0{i} * Hy(seg,:);       % solution
        Sigma_x{i} = Sigma0{i} - Sigma0{i} * HPhi(seg,seg) * Sigma0{i};
         
        
        currentLoc = seg(end);
              
        %=========== Learn intra-frame structure in blocks ===========
        if LEARNTYPE == 0
            A{i} = eye(currentLen);
            invA{i} = eye(currentLen);
        % regularization
        elseif LEARNTYPE == 1
            if equalSize == 0
                if currentLen > 1
                    temp = Sigma_x{i}/gamma(i) + mu_x(seg,:)*mu_x(seg,:)'/L/gamma(i);
                    
                    r0 = r0 + mean(diag(temp));
                    r1 = r1 + mean(diag(temp,1));
                end
            elseif equalSize == 1
                A0 = A0 + Sigma_x{i}/gamma(i)/L + mu_x(seg,:)*mu_x(seg,:)'/L/gamma(i);
            end

        end % end of learnType

    end % end of usedNum
    
    %=========== Learn inter-frame structures with Constraint ===========
    % If blocks have the same size
    if (equalSize == 1) & (LEARNTYPE == 1)

        % Toeplitz structure to regularize the covariance structure
        b = (mean(diag(A0,1))/mean(diag(A0)));
        if abs(b) >= 0.90, b = 0.90*sign(b); end;
        bs = [];
        for j = 1 : maxLen, bs(j) = (b)^(j-1); end;
        A0 = toeplitz(bs);
 
        for i = 1 : usedNum
             
            A{i} = A0;
            invA{i} = inv(A{i});
        end
    
    % if blocks have different sizes
    elseif (equalSize == 0) & (LEARNTYPE == 1)
        r = r1/r0; if abs(r) >= 0.99, r = 0.99*sign(r); end;

        for i = 1 : usedNum
            currentLen = size(Sigma_x{i},1);

            bs = [];
            for j = 1 : currentLen, bs(j) = r^(j-1); end;
            A{i} = toeplitz(bs);
            invA{i} = inv(A{i});
        end 
    end

    
    % estimate gamma(i) and lambda 
    if LearnLambda == 1         % this learning rule is suitable for general noisy cases (SNR<=20dB)
        gamma_old = gamma;
        lambdaComp = 0; 
        currentLoc = 0;
        invSigma_y = eye(M)/(lambda * eye(M) + PhiAPhi);
        
        for i =  1 : usedNum

            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            lambdaComp = lambdaComp + trace(Phi(:,currentSeg)*Sigma_x{i}*Phi(:,currentSeg)');
            currentLoc = currentSeg(end);

            g_num = sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/L;
            g_den = trace(invSigma_y * Phi(:,currentSeg) * A{i} * Phi(:,currentSeg)');
            gamma(i) = sqrt(g_num/g_den);
            
            Sigma0{i} = A{i} * gamma(i);
            
        end
        lambda = norm(y - Phi*mu_x,'fro')^2/(M*L) + lambdaComp/M;
        
    elseif LearnLambda == 2   % this learning rule is suitable for high SNR cases (>20dB)
        gamma_old = gamma;
        lambdaComp = 0;  
        currentLoc = 0;
        invSigma_y = eye(M)/(lambda * eye(M) + PhiAPhi);
        for i =  1 : usedNum                 
            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            
            lambdaComp = lambdaComp + trace(Sigma_x{i}* Phi(:,currentSeg)'*Phi(:,currentSeg));
            %====================original method===============================
%             gamma(i) = trace(invA{i}*Sigma_x{i})/currentLen ...
%                 + sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/(L*currentLen);
            %=============== Bound-optimization method==========================            
            g_num = sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/L;
            g_den = trace(invSigma_y * Phi(:,currentSeg) * A{i} * Phi(:,currentSeg)');
            gamma(i) = sqrt(g_num/g_den);
            
            Sigma0{i} = A{i} * gamma(i);
            
            currentLoc = currentSeg(end);
        end
%         lambda = norm(y - Phi * mu_x,'fro')^2/(M*L) + lambda * (length(mu_x) - lambdaComp)/M; 
        lambda = norm(y - Phi * mu_x,'fro')^2/(M*L) + lambdaComp/M; 
    else   % only estimate gamma(i)
        gamma_old = gamma;
        currentLoc = 0;
        invSigma_y = eye(M)/(lambda * eye(M) + PhiAPhi);
        for i =  1 : usedNum
            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            
            g_num = sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/L;
            g_den = trace(invSigma_y * Phi(:,currentSeg) * A{i} * Phi(:,currentSeg)');
            gamma(i) = sqrt(g_num/g_den);
            
            Sigma0{i} = A{i} * gamma(i);
            
            currentLoc = currentSeg(end);
        end
    end

    
    % ================= update B ============
    mu_old = mu_t;
    mu_t = mu_x * sqrtm(B); 
    B = zeros(L); currentLoc = 0;
    for i = 1 : usedNum
      
        currentLoc = currentLoc + 1;
        currentLen = size(Sigma0{i},1);
        seg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        B = B + mu_t(seg,:)' * invA{i} * mu_t(seg,:)/gamma_old(i) + MatrixReg;  
        
        currentLoc = seg(end);
    end
    B = B/norm(B,'fro');
    if size(Phi,2) <= size(Phi,1), B = eye(L); end;
     
    % ================= Check stopping conditions ==============
        % check convergence
        dmu_change = dmu; %%%%% add extra condition
        if (size(mu_t) == size(mu_old))
            dmu = max(max(abs(mu_old - mu_t)));
            if (dmu < EPSILON)  break;  end;
        end;
        if (dmu==dmu_change) count_num= count_num+1; end %%% add 
        if (PRINT)
            disp([' iters: ',num2str(count),...
                ' num coeffs: ',num2str(usedNum), ...
                ' min gamma: ', num2str(min(gamma)),...
                ' gamma change: ',num2str(max(abs(gamma - gamma_old))),...
                ' mu change: ', num2str(dmu)]);
        end;
        if (count >= MAX_ITERS), if PRINT, fprintf('Reach max iterations. Stop\n\n'); end; break;  end;
        if (count_num >= 6), if PRINT, fprintf('No change. Stop\n\n'); end; break; end  %%%%% 
    
end;
%% Expand hyperparameyers
gamma_used = sort(keep_list);
gamma_est = zeros(p,1);
gamma_est(keep_list,1) = gamma;


%% reconstruct the original signal
x = zeros(N,L);
currentLoc = 0;
for i = 1 : usedNum

    currentLen = size(Sigma0{i},1);
    currentLoc = currentLoc + 1;
    seg = currentLoc : 1 : currentLoc + currentLen - 1;

    realLocs = blkStartLoc0(keep_list(i)) : blkStartLoc0(keep_list(i))+currentLen-1;

    x( realLocs,: ) = mu_t( seg,: );
    currentLoc = seg(end);
end

if (scl < 0.4) | (scl > 1)
    Result.x = x * scl/0.4;
else
    Result.x = x;
end

Result.gamma_used = gamma_used;
Result.gamma_est = gamma_est;
Result.A = A;
Result.count = count;
Result.lambda = lambda;
return;

