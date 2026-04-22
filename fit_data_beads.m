%% Load the raw data
% notes=importdata('N:\share\CHT\TrpASensor\trpasensornote_all.txt');
% notes=importdata('N:\share\CHT\TrpASensor\ctrlsensornote.txt');
% notes=importdata('N:\share\CHT\TrpASensor\TELC_notes.txt');
addpath(genpath('W:\AccuSleep\BEADS_toolbox'))
notes=importdata('N:\share\CHT\TrpASensor\TELC_notes.txt');
% Filter parameters

tic

% fc : cut-off frequency (cycles/sample)

d = 1;          % d : filter order parameter (d = 1 or 2)
% Positivity bias (peaks are positive)
r = 12;          % r : asymmetry parameter
% Regularization parameters
amp = 0.00015; % 0.006
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

% start=0; % h
% all_notes_mean_std=nan(length(notes),2);
% dF_labels_all=nan(length(notes),3);
% dF_labels_norm_all=nan(length(notes),3);
% dF_labels_beads_all=nan(length(notes),3);
% dF_labels_beads_norm_all=nan(length(notes),3);

load('W:\newdata_for_revision\fiberphotometry_newdata\TrpA_sensor_LH\s333-250508-145646\1\info.mat')
tic
% list=1:length(notes);
% list([63 72 82 93])=[];
%  l=1:length(fc_all)

% list=dir();

fc = 0.09;%fc_all(l);
for j=1:length(notes)

%     temp=list(j).name;

%     if isfolder(temp)
%         for m=1:2
%             path = strcat(list(j).folder,'\',list(j).name,'\',num2str(m));


                path = notes{j};
%             if isfile(strcat(path,'\dF_labels.mat'))
%                 continue
%             end
            try
                load (strcat(path,'\labels.mat'));
                load (strcat(path,'\gcamp.mat'));
%                 load (strcat(path,'\info.mat'));
                %         load (strcat(path,'\dF_labels.mat'));
            catch
                continue
            end
            %     dF_labels_beads_009=dF_labels_beads;
            %
            dF=double(dF);
            start=0;
            stop=length(dF)/info.samplerate2/3600;

            dF=dF(1+round(start*3600*info.samplerate2):end);
            labels=labels(start*1440+1:min(end,round(stop*1440)));

            idx_labels=round((0:length(labels)-1).*info.samplerate2.*2.5);
            idx_labels_full=nan(length(idx_labels),2543); % 2.5s sr
            for i=1:length(idx_labels)
                idx_labels_full_temp=dF(idx_labels(i)+1:min(idx_labels(i)+2543,end));
                idx_labels_full(i,1:length(idx_labels_full_temp))=idx_labels_full_temp;
            end
            dF_labels=nanmedian(idx_labels_full,2);
            dF_labels2=dF_labels*-1;

            [x1, f1, cost] = beads(dF_labels2, d, fc, r, lam0, lam1, lam2);

            dF_labels_bl=f1*-1;
            dF_labels_beads=x1*-1;


            idx_labels=round((0:(length(labels)-1)*5).*info.samplerate2.*0.5);
            idx_labels_full=nan(length(idx_labels),508); %0.5s sr
            for i=1:length(idx_labels)
                idx_labels_full_temp=dF(idx_labels(i)+1:min(idx_labels(i)+508,end));
                idx_labels_full(i,1:length(idx_labels_full_temp))=idx_labels_full_temp;
            end
            dF_labels_05=nanmedian(idx_labels_full,2);

            dF_labels2_05=dF_labels_05*-1;

            [x1, f1, cost] = beads(dF_labels2_05, d, fc, r, lam0, lam1, lam2);
            dF_labels_05_bl=f1*-1;
            dF_labels_05_beads=x1*-1;

            save(strcat(path,'\dF_labels.mat'),'dF_labels','dF_labels_05','dF_labels_bl','dF_labels_beads')
            save(strcat(path,'\dF_labels.mat'),'dF_labels_beads','dF_labels_05_bl','dF_labels_05_beads',"-append")


            % time=1:length(dF_labels);
            % time=time./1440;
            %  figure
            %  plot(time,dF_labels,'LineWidth',1)
            %  hold on
            %  plot(time,dF_labels_beads,'LineWidth',1)
            %
            %  y_max=max([dF_labels_beads;dF_labels]);
            %  y_min=min([dF_labels_beads;dF_labels]);
            %  labels1=labels;
            % %  labels1=repmat(labels',5,1);
            % %  labels1=labels1(1:end);
            %  for i=1:length(labels1)
            %          v = [(i-1)*info.epochtime/3600 y_min;(i-1)*info.epochtime/3600+info.epochtime/3600 y_min;(i-1)*info.epochtime/3600+info.epochtime/3600 y_max;(i-1)*info.epochtime/3600 y_max];
            %          f = [1 2 3 4];
            %         if labels1(i)==1
            %             patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceColor',[0.9686 0.5765 0.1176],'FaceAlpha',0.3);
            %         elseif labels1(i)==2
            %             patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
            %             elseif labels1(i)==3
            %             patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceColor',[0.4 0.7608 0.5608],'FaceAlpha',0.3);
%         end
%     end



end
toc
%%
function [x, f, cost] = beads(y, d, fc, r, lam0, lam1, lam2)

% [x, f, cost] = beads(y, d, fc, r, lam0, lam1, lam2)
%
% Baseline estimation and denoising using sparsity (BEADS)
%
% INPUT
%   y: Noisy observation
%   d: Filter order (d = 1 or 2)
%   fc: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5)
%   r: Asymmetry ratio
%   lam0, lam1, lam2: Regularization parameters
%
% OUTPUT
%   x: Estimated sparse-derivative signal
%   f: Estimated baseline
%   cost: Cost function history

% Reference:
% Chromatogram baseline estimation and denoising using sparsity (BEADS)
% Xiaoran Ning, Ivan W. Selesnick, Laurent Duval
% Chemometrics and Intelligent Laboratory Systems (2014)
% doi: 10.1016/j.chemolab.2014.09.014
% Available online 30 September 2014

% The following parameter may be altered.
Nit = 30;       % Nit: Number of iterations
pen = 'L1_v2';  % pen : penalty function for sparse derivative ('L1_v1' or 'L1_v2')
EPS0 = 1e-6;    % cost smoothing parameter for x (small positive value)
EPS1 = 1e-6;    % cost smoothing parameter for derivatives (small positive value)

switch pen
    case 'L1_v1'
        phi = @(x) sqrt(abs(x).^2 + EPS1);
        wfun = @(x) 1./(sqrt(abs(x).^2 + EPS1));
    case 'L1_v2'
        phi = @(x) abs(x) - EPS1 * log(abs(x) + EPS1);
        wfun = @(x) 1./( abs(x) + EPS1);
    otherwise
        disp('penalty must be L1_v1, L1_v2')
        x = []; cost = []; f = [];
        return
end

theta = @(x) sum(x(x>EPS0)) - r * sum(x(x<-EPS0)) ...
    + sum( (1+r)/(4*EPS0)*x(abs(x)<=EPS0).^2 ...
    + (1-r)/2 * x(abs(x)<=EPS0) + EPS0*(1+r)/4 );

y = y(:);
x = y;
cost = zeros(1, Nit);
N = length(y);
[A, B] = BAfilt(d, fc, N);
H = @(x) B*(A\x);
e = ones(N-1, 1);
D1 = spdiags([-e e], [0 1], N-1, N);
D2 = spdiags([e -2*e e], 0:2, N-2, N);
D = [D1;  D2];
BTB = B'*B;

w = [lam1 * ones(N-1, 1); lam2 * ones(N-2, 1)];
b = (1-r)/2 * ones(N, 1);
d = BTB * (A\y) - lam0 * A' * b;

gamma = ones(N, 1);

for i = 1:Nit

    Lambda = spdiags(w.*wfun(D*x), 0, 2*N-3, 2*N-3);

    k = abs(x) > EPS0;
    gamma(~k) = ((1 + r)/4) / abs(EPS0);
    gamma(k) = ((1 + r)/4) ./  abs(x(k));
    Gamma = spdiags(gamma, 0, N, N);

    M = 2 * lam0 * Gamma + D' * Lambda * D;
    x = A * ((BTB + A'*M*A)\d);

    cost(i) = 0.5 * sum(abs(H(y - x)).^2) + lam0 * theta(x) ...
        + lam1 * sum(phi(diff(x))) + lam2 * sum(phi(diff(x, 2)));
end

f = y - x - H(y-x);

end


% --- local function ----

function [A, B] = BAfilt(d, fc, N)
% [A, B] = BAfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass filter.
% The matrices are 'sparse' data type in MATLAB.
%
% INPUT
%   d  : degree of filter is 2d (use d = 1 or 2)
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal

b1 = [1 -1];
for i = 1:d-1
    b1 = conv(b1, [-1 2 -1]);
end
b = conv(b1, [-1 1]);

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;

a = 1;
for i = 1:d
    a = conv(a,[1 2 1]);
end
a = b + t*a;

A = spdiags( a(ones(N, 1), :), -d:d, N, N);   % A: Symmetric banded matrix
B = spdiags(b(ones(N, 1), :), -d:d, N, N);    % B: banded matrix

end

