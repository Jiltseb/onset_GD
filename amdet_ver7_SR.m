%Copyright (C) 2017 Speech and Music Technology Lab,
%Indian Institute of Technology Madras
                
%This file is part of GD based onset detection(onset_GD).                         
%onset_GD is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or   
%(at your option) any later version.
                
%This software is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.                             

%You should have received a copy of the GNU General Public License          
%.  If not, see <http://www.gnu.org/licenses/>

% AM Demoduation and Grp delay based onset detection 
% AUTHOR :  P A Manoj Kumar
% Edited by : Jilt Sebastian
% Note: The placement of the file is essential. The binary WordSegmentWithSielnceRemoval and the script extrema.m must be in the
% same folder. This is the single resolution program, hence wsF is a parameter!



%clear all;
close all;
%clc;
Mlf_file = 'result_1.mlf';

% matlabpool(12);
warning('off','all');
warning;
% Parameter definitions

downsampling_rate = 10:2:10;
smoothening_factor = 44:2:44;    
winScaleFactor = 30:5:30;

thres = 0.2:0.02:0.2;

F_measure = zeros(length(downsampling_rate),length(smoothening_factor),length(thres));


cd test_here;
listing = dir(pwd);
cd ..;


% Begin of paramterization
for indexA = 1:length(downsampling_rate)
    for indexB = 1:length(smoothening_factor)
            for indexD = 1:length(thres)
                fprintf('Paramterization details : dwnsmpl_rate = %d\n smth_factor = %d, thres = %f\n',downsampling_rate(indexA),smoothening_factor(indexB),thres(indexD));
%                     fid2 = fopen(Number_alone_file,'w');
                fid3 = fopen(Mlf_file,'w');
                fprintf(fid3,'#!MLF!#\n');
                for ii = 3:1:length(listing)

                    %======================================================================
                    % Part1: Using Amplitude Demodulation, and applying Group delay on it
                    filename = listing(ii).name;
                    orig_filename = filename;
                    cd test_here;
                    [Y,Fs]=wavread(filename);
                    cd ..;    

                    DF=diff(Y);                 % Differentiate it to emphasize the frequency components in amplitude
                    fprintf('Working on %s\n',filename);

                    a = hilbert(DF);
                    Z=abs(a+DF);

                    D = downsample(Z,downsampling_rate(indexA));
                    S = smooth(D,smoothening_factor(indexB),'moving');   % This is the measure!
                    assignin('base','S',S);
                    assignin('base','Y',Y);
                    grp_delay = ones(length(S),1);
                    gd_sum = ones(length(S),1);

                    parfor wsfIndex = 1:length(winScaleFactor)
                        
                        tempDir = sprintf('temp_%d',wsfIndex);
                        mkdir(tempDir); cd(tempDir);
                        energy_file_name = strcat(filename(1:end-4),'.en');
                        dlmwrite(energy_file_name,S*1000,'\n');
                        spec_file_name = energy_file_name(1:end-2);
                        spec_file_name =strcat(spec_file_name,'spec');
                        % Invoking the binary                        
                        copyfile('../fe-words.base_ref','fe-words.base');
                        ctrl_file = 'fe-words.base';
                        temp_ctrl_file = strcat('temp.base');
                        % Changing the winscalefactor parameter in config file
                        a = importdata(ctrl_file);
                        a = struct2cell(a);
                        a{1}(3) = winScaleFactor(wsfIndex);
                        fprintf('Window scale factor is %d\n',winScaleFactor(wsfIndex));
                        fid0 = fopen(temp_ctrl_file,'w');
                        for i = 1:length(a{1})
                            fprintf(fid0,'%s %s %f\n',char(a{2}(i,1)),char(a{2}(i,2)),a{1}(i));
                        end
                        copyfile(temp_ctrl_file,ctrl_file);
                        delete(temp_ctrl_file);
                        fclose(fid0);

                        dummy1 = 'b';
                        dummy2 = 'c';
                        dummy3 = 'd';
                        dummy4 = 'e';
                        dump = 'dump.txt';

                        system(sprintf('../WordSegmentWithSilenceRemoval %s %s %s %s %s %s %s > %s 2>&1',ctrl_file,energy_file_name,spec_file_name,dummy1,dummy2,dummy3,dummy4,dump));

                        delete(energy_file_name);
                        temp = load(spec_file_name);
                        delete(spec_file_name);
                        temp = temp(:,5);
                        temp(length(S)+1:end) = [];
%                         temp = smooth(temp,2*smoothening_factor(indexB),'moving');
                        grp_delay = grp_delay.*temp;
                        temp = temp - mean(temp);
                        gd_sum = gd_sum + cumsum(temp);
%                         plot(temp);
%                         eval(sprintf('temp_%d = temp;',j));
%                         figure;
                        cd ..; %rmdir(tempDir,'s');
                    end
%                     gd_sum = smooth(gd_sum,2*smoothening_factor(indexB),'moving');
                    grp_delay = diff(gd_sum);
                    grp_delay = smooth(grp_delay,2*smoothening_factor(indexB),'moving');   % A moving average with 1 ms interval
                    
                    grp_delay = grp_delay/max(grp_delay);                                        
                    assignin('base','grp_delay',grp_delay);
                    %======================================================================



                    %======================================================================
                    % Part2: Reading the contents of group delay file, and getting the
                    % onsets
                    threshold = thres(indexD);

                    stroke_loc = zeros(1,length(grp_delay));
                    % Go to each minima, and calculate height till next maxima. Keep a
                    % threshold on this to decide if stroke!


                    t = 1:length(grp_delay);
            %         figure
            %         plot(t,grp_delay)
                    [ymax,imax,ymin,imin] = extrema(grp_delay);
            %         hold on
            %         plot(t(imax),ymax,'r*',t(imin),ymin,'g*')

                    % sort the minimas and maximas;
                    temp_min = sortrows([imin ymin]);
                    imin = temp_min(:,1)';
                    ymin = temp_min(:,2)';
                    clear temp_min;

                    temp_max = sortrows([imax ymax]);
                    imax = temp_max(:,1)';
                    ymax = temp_max(:,2)';
                    clear temp_max;



                    if (imin(1) < imax(1) )  % fine, just truncate the maximum
                        imin(1) = []; ymin(1) = [];

                        if (length(imin) > length(imax) )
                            imin(length(imax)+1:end) = [];
                            ymin(length(imax)+1:end) = [];
                        elseif (length(imin) < length(imax) )
                            imax(length(imin)+1:end) = [];
                            ymax(length(imin)+1:end) = [];
                        end
                    else                                                    

%                             imax(1) = []; ymax(1) = [];
                        if (length(imin) > length(imax) )
                            disp('this shouldnt have come');
                            imin(length(imax)+1:end) = [];
                            ymin(length(imax)+1:end) = [];
                        elseif (length(imin) < length(imax) )

                            imax(length(imin)+1:end) = [];
                            ymax(length(imin)+1:end) = [];
                        end
                    end


            %         if (imin(1) == 1)
            %             imin(1) = []; ymin(1) = []; % Making them of equal size
            %             imax(end) = []; ymax(end) = [];
            %         end
            %         if (imax(1) == 1)
            %             imax(1) = []; ymax(1) = [];
            %             imin(end) = []; ymin(end) = [];
            %         end

                    assignin('base','ymax',ymax);
                    assignin('base','imax',imax);
                    assignin('base','ymin',ymin);
                    assignin('base','imin',imin);
                    assignin('base','grp_delay',grp_delay);

%                     fprintf('the differntiating parameter for file %s is %f\n',filename,-mean(ymax(ymax>0))/mean(ymin(ymin<0)));

                     %==================================================================
                     % Algorithm1  for stroke location
                     index_stroke = 1;
                     peak_valley_heights = ymax - ymin;

                     for index = 1:1:length(peak_valley_heights)
                         if (peak_valley_heights(index) > threshold)
                     %         fprintf('Interest point at %d\n',ceil((imin(index) + imax(index))/2)); 
                             stroke_loc(index_stroke) = ceil((imin(index) + imax(index))/2);
                             index_stroke = index_stroke + 1;
                         end
                     end
                     %==================================================================


                    stroke_loc(stroke_loc==0) = [];

                    assignin('base','stroke_loc',stroke_loc);
                    assignin('base','peaks',peaks);

                    %======================================================================



                    %======================================================================
                    % Printing in standard MLF format
                    dangerflag = 0;
                    cd test_here;
                    [ X, Fs] = wavread(orig_filename);
                    cd ..;

%                         fid2 = fopen(Number_alone_file,'a');
                    fid3 = fopen(Mlf_file,'a');

                    length_wav_file = length(X)*1/Fs;
                    stroke_loc = stroke_loc*downsampling_rate(indexA)/Fs;  % Converting into seconds
                    if (isempty(stroke_loc))        % Provision for null strokes
                        filename = filename(1:end-4); 
                        fprintf(fid3,'\"*/%s.lab\"\n0\t%f\n%f\t%f',filename,length_wav_file-0.0100000,length_wav_file-0.0100000,length_wav_file);
                    else                   
                        if (length_wav_file - 0.01 < stroke_loc(end) )  % Artificial stroke at end due to compuation of group delay function
                            stroke_loc(end) = [];
                            if (isempty(stroke_loc))        % Provision for only one stroke, that too at end of file
                                filename = filename(1:end-4); 
                                fprintf(fid3,'\"*/%s.lab\"\n0\t%f\n%f\t%f',filename,length_wav_file-0.0100000,length_wav_file-0.0100000,length_wav_file);
                                dangerflag = 1;
                            end
                        end
                        if (dangerflag~=1)
                            clear X;
                            filename = filename(1:end-4);
                            assignin('base','filename',filename);
                            assignin('base','stroke_loc',stroke_loc);
                            fprintf(fid3,'\"*/%s.lab\"\n0\t%f\n',filename,stroke_loc(1));
                            index3 = 1;
                            while (index3 < length(stroke_loc) )
                                fprintf(fid3,'%f\t%f\n',stroke_loc(index3),stroke_loc(index3 + 1));
                                index3 = index3 + 1;
                            end
                            fprintf(fid3,'%f\t%f\n',stroke_loc(index3),length_wav_file-0.0100000);  % Forcing 10ms for silence at the end
                            fprintf(fid3,'%f\t%f\n',length_wav_file-0.0100000,length_wav_file);  
                        end
                    end

%                         fprintf(fid2,'Number of strokes in file %s is %d\n',listing(ii).name, length(stroke_loc));
%                         fclose(fid2);
                    fclose(fid3);

                    %====================================================================== 
                end
                
            end
    end
end
