function [dwp_mat, dwm_mat, other] = errorPosition_20240606(U, Xref, K, Wp, Wm, M, tracking_error)
% Prep
N = size(U,1);
tau = size(U,3);
Imag = 0+1j;
Delta = @(a,b) double(logical(a==b));
% Approximate du and dxref from tracking error
dxref_mat = ones(size(Xref)).*tracking_error./sqrt(tau);
du_mat = ones(size(Xref)).*tracking_error;

% Compute error, with matrices
D_mat = zeros(size(K,1),2,2); % K*2*2
Tr_mat = zeros(size(K,1),1); % K*1
Det_mat = zeros(size(K,1),1); % K*1
wp_mat = Wp; % K*1, \omega_plus, calculated
wm_mat = Wm; % K*1, \omega_minus, calculated

dwp_dTr_mat = zeros(size(K,1),1); % K*1
dwm_dTr_mat = zeros(size(K,1),1); % K*1
dwp_dDet_mat = zeros(size(K,1),1); % K*1
dwm_dDet_mat = zeros(size(K,1),1); % K*1

dwp_dxref_mat = zeros(size(K,1), N, 2); % K*N*2
dwm_dxref_mat = zeros(size(K,1), N, 2); % K*N*2
dwp_du_mat = zeros(size(K,1), N, 2, tau); % K*N*2*T
dwm_du_mat = zeros(size(K,1), N, 2, tau); % K*N*2*T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dD_du_mat = zeros(N, 2,  size(K,1), 2, 2, tau); % prev N*2*T*K*2*2 now N*2*K*2*2*T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dD_du_mat = zeros(N, 2, tau, size(K,1), 2, 2); % N*2*T*K*2*2
dD_dxref_mat = zeros(N, 2, size(K,1), 2, 2); % N*2*K*2*2

% dxref_mat = zeros(N, 2); % N*2
% du_mat = zeros(N,2); % N*2, [the same to any t, so no such axis]

% dwp_mat = zeros(size(K,1),1); % K*1
% dwm_mat = zeros(size(K,1),1); % K*1

% Get D-1 matrix 'D_mat'
for i = 1:size(K,1)
    for j = 1:2
        for k = 1:2
            D_mat(i,j,k) = M(K(i,1),K(i,2), j, k);
        end; end; end

% Get Trace and Detterminant 'Tr_mat', 'Det_mat'
for i = 1:size(K,1)
    Tr_mat(i) = D_mat(i,1,1) + D_mat(i,2,2);
    Det_mat(i) = D_mat(i,1,1)*D_mat(i,2,2) - D_mat(i,1,2)*D_mat(i,2,1);
end

% Get dw_dTr and dw_dDet
for i = 1:size(K,1)
    dwp_dTr_mat(i) = -(1/4) * wp_mat(i)^3 * (1 - Tr_mat(i) / sqrt( (Tr_mat(i)/2)^2 - Det_mat(i) ) );
    dwm_dTr_mat(i) = -(1/4) * wm_mat(i)^3 * (1 + Tr_mat(i) / sqrt( (Tr_mat(i)/2)^2 - Det_mat(i) ) );
    dwp_dDet_mat(i) = -(1/4) * wp_mat(i)^3 * ( + 1 / sqrt( (Tr_mat(i)/2).^2 - Det_mat(i) ) );
    dwm_dDet_mat(i) = -(1/4) * wm_mat(i)^3 * ( - 1 / sqrt( (Tr_mat(i)/2).^2 - Det_mat(i) ) );
end

% Get dD_dxref and dD_du
figure(1);clf;
timing_take=tic;
hold on
loopnumb = 0;
for i = 1:N % k
    if ~mod(i,10)
        disp(['Error estimation: ', num2str(i), ' / ', num2str(N)])
        
    end
    dXref_temp = Xref(i,:) - Xref;
    for ii = 1:2 % a
        for jj = 1:2 % b
            u_ia_kb_t = reshape(U(:,ii,:),N,tau) * reshape(U(i,jj,:),tau,1);
            u_ib_ka_t = reshape(U(:,jj,:),N,tau) * reshape(U(i,ii,:),tau,1);
            for j = 1:2 % g
                delta_ag = Delta(j,ii);
                delta_bg = Delta(j,jj);
                for l = 1:size(K,1)
                    
                    qgamma_temp = K(l,j);
                    dD_dxref_mat(i,j,l,ii,jj) = Imag * qgamma_temp / N / tau *...
                        ( sum( exp(Imag * dXref_temp * K(l,:)') .* u_ia_kb_t ) -...
                        sum( exp(Imag * -dXref_temp * K(l,:)') .* u_ib_ka_t ));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start
                    temp_k = zeros([1 tau]);
                    for k = 1:tau % t
                        firsttemp = exp(Imag * dXref_temp * K(l,:)');
                        firstsum = firsttemp'*U(:,ii,k);
                        
                        secondtemp = exp(Imag * -dXref_temp * K(l,:)');
                        secondsum = secondtemp'*U(:,jj,k);
                        %                         dD_du_mat(i,j,k,l,ii,jj) = 1 / N / tau *...
                        %                             ( delta_bg * firstsum  + ...
                        %                             delta_ag * secondsum);
                        temp_k(k) = delta_bg*firstsum+ delta_ag*secondsum;
                        loopnumb = loopnumb + 1;
                    end
                    dD_du_mat(i,j,l,ii,jj,:) = 1 / N / tau *temp_k;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end
                    
                end
                
%                 %%%%%%%%%for plotting the expected runtime%%%%%%%%%%%%%%
%                 disp([num2str(loopnumb) '/' num2str(N*2*2*2*size(K,1)*tau)]);
%                 timing_end=toc(timing_take);
%                 scatter(loopnumb,timing_end/loopnumb*(N*2*2*2*size(K,1)*tau),30,'k','filled')
%                 scatter(loopnumb,timing_end,30,'r','filled')
%                 drawnow
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end; end; end
    
end


% Get dw_dxref and dw_du
for i = 1:size(K,1) % k1,k2
    for j = 1:N % i
        for k = 1:2 % a
            part1_temp = (dD_dxref_mat(j,k,i,1,1) + dD_dxref_mat(j,k,i,2,2));
            part2_temp = (D_mat(i,2,2) * dD_dxref_mat(j,k,i,1,1) + ...
                D_mat(i,1,1) * dD_dxref_mat(j,k,i,2,2) - ...
                D_mat(i,1,2) * dD_dxref_mat(j,k,i,2,1) - ...
                D_mat(i,2,1) * dD_dxref_mat(j,k,i,1,2));
            dwp_dxref_mat(i,j,k) = dwp_dTr_mat(i) * part1_temp + ...
                dwp_dDet_mat(i) * part2_temp;
            dwm_dxref_mat(i,j,k) = dwm_dTr_mat(i) * part1_temp + ...
                dwm_dDet_mat(i) * part2_temp;
            
            for l = 1:tau
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start
                part1_t_temp = (dD_du_mat(j,k,i,1,1,l) + dD_du_mat(j,k,i,2,2,l));
                part2_t_temp = (D_mat(i,2,2) * dD_du_mat(j,k,i,1,1,l) + ...
                    D_mat(i,1,1) * dD_du_mat(j,k,i,2,2,l) - ...
                    D_mat(i,1,2) * dD_du_mat(j,k,i,2,1,l) - ...
                    D_mat(i,2,1) * dD_du_mat(j,k,i,1,2,l));
                dwp_du_mat(i,j,k,l) = dwp_dTr_mat(i) * part1_t_temp +...
                    dwp_dDet_mat(i) * part2_t_temp;
                dwm_du_mat(i,j,k,l) = dwm_dTr_mat(i) * part1_t_temp +...
                    dwm_dDet_mat(i) * part2_t_temp;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end
            end
        end; end; end

% Get final dw
dxref_mat_repeat = reshape(dxref_mat, 1,N,2);
dxref_mat_repeat = repmat(dxref_mat_repeat, size(K,1),1,1);
dwp_part1_temp = (dwp_dxref_mat.^2) .* (dxref_mat_repeat.^2);
dwm_part1_temp = (dwm_dxref_mat.^2) .* (dxref_mat_repeat.^2);

du_mat_repeat = reshape(du_mat, 1,N,2,1);
du_mat_repeat = repmat(du_mat_repeat, size(K,1),1,1,tau);
dwp_part2_temp = (dwp_du_mat.^2) .* (du_mat_repeat.^2);
dwm_part2_temp = (dwm_du_mat.^2) .* (du_mat_repeat.^2);

dwp_mat = sqrt(sum(dwp_part1_temp, [2,3]) + sum(dwp_part2_temp, [2,3,4]));
dwm_mat = sqrt(sum(dwm_part1_temp, [2,3]) + sum(dwm_part2_temp, [2,3,4]));

dwp_mat_part1 = sqrt( sum(dwp_part1_temp, [2,3]));
dwm_mat_part1 = sqrt( sum(dwm_part1_temp, [2,3]));
dwp_mat_part2 = sqrt( sum(dwp_part2_temp, [2,3,4]));
dwm_mat_part2 = sqrt( sum(dwm_part2_temp, [2,3,4]));

other.dwp_mat_part1 = dwp_mat_part1;
other.dwp_mat_part2 = dwp_mat_part2;
other.dwm_mat_part1 = dwm_mat_part1;
other.dwm_mat_part2 = dwm_mat_part2;
end