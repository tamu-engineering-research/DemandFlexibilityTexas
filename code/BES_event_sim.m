%% intro
clear;
define_constants;
case_path = "..\case\case_Mar3_5pm.mat";
gen_path = "..\..\blackout\From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)\ERCOT Generation by Source\2021-02-01 to 2021-02-19.csv";
load_path = "..\..\blackout\From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)\ERCOT Demand by Subregion\2021-02-01 to 2021-02-19.csv";
wind_path = "..\..\blackout\From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)\BES profiles\BES_wind_Feb2021_v20210224.csv";
solar_path = "..\..\blackout\From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)\BES profiles\BES_solar_Feb2016_v20210224.csv";
hydro_path = "..\..\blackout\From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)\BES profiles\BES_hydro_Feb2016_v20210224.csv";
lf_path =  "..\..\blackout\From EIA (Generation, Demand, Demand Forecasts, and Interchange by BA)\ERCOT Demand, Demand Forecast, and Net Interchange\2021-02-01 to 2020-02-19.csv";
out_ng_path = "..\..\blackout\data\unit_outage\ng_outage.csv";
out_pv_path = "..\..\blackout\data\unit_outage\pv_outage.csv";
out_coal_path = "..\..\blackout\data\unit_outage\coal_outage.csv";
gen_county_path = "..\..\blackout\data\unit_outage\gen_county.csv";
ng_derate_path = "..\..\blackout\data\unit_outage\ng_derate_fac.csv";

out_path = "..\..\blackout\data\out_mat.mat";
shed_path = "..\..\blackout\data\shed_mat.mat";
cap_path = "..\..\blackout\data\ercotcap_mat.mat";
lcf_path = "..\..\blackout\data\lf_mat.mat";
cty_outage_pct_path = '..\..\blackout\data\load_outage\cty_out_pct.csv';
cty_outage_dist_path = '..\..\blackout\data\load_outage\cty_out_dist.csv';
all_bus_cty_path = '..\..\blackout\data\load_outage\bus_county.csv';

% ERCOT load resources
ers_path = '..\data\ERS(Feb-12-18).csv';
rrs_path = '..\data\Load Resourses (Feb-12-18).csv';

[ers_ts, ~, ~] = xlsread(ers_path);
[rrs_ts, ~, ~] = xlsread(rrs_path);


% load composition and profiling
bus_pop_sec_path = '..\data\bus_pop_sec.csv';
area_profile_folder = '..\data\area_profiles_pct\';

[bus_pop, ~, ~] = xlsread(bus_pop_sec_path);


% load profiling
area_mpc = ["FW", "N", "W", "S", "NC", "SC", "C", "E"];
area_profiles_pct = containers.Map;

for key = area_mpc
    area_profile_path = append(area_profile_folder, key, ".csv");
    [area_profiles_pct(key), ~, ~] = xlsread(area_profile_path);
end


out_ts = load(out_path);
out_ts = out_ts.out_mat;
shed_ts = load(shed_path);
shed_ts = shed_ts.shed_mat;
shed_ts(1:312,3) = NaN;
cap_ts = load(cap_path);
cap_ts = cap_ts.cap_mat;
lcf_ts = load(lcf_path);
lcf_ts = lcf_ts.lf_mat;

mpc = loadcase(case_path);
bus_num = size(mpc.bus,1);
all_bus_pd = mpc.bus(:, PD);

% configs
mpopt = mpoption('pf.nr.max_it', 50);
mpopt = mpoption(mpopt,'out.all',0);
mpopt = mpoption(mpopt,'verbose',0); 

[gen_val, gen_ts, ~] = xlsread(gen_path);
gen_ts(1, :) = [];
gen_ts = gen_ts(:, 2);
gen_val(1:24,:) = [];
gen_val(end,:) = [];

[load_val, load_ts, ~] = xlsread(load_path);
load_ts(1, :) = [];
load_val(1:24, :) = [];

[lf_val, lf_ts, ~] = xlsread(lf_path);
lf_ts(1, :) = [];
lf_val(1:24, :) = [];

% load out data from 3rd party
[~, bus_county, ~] = xlsread(all_bus_cty_path);

[cty_out_pct_val, cty_out_pct_ts, ~] = xlsread(cty_outage_pct_path);
cty_out_header = cty_out_pct_ts(1,2:end);
bus_county_num = size(cty_out_header,2);
cty_out_pct = [zeros(11*24, bus_county_num); cty_out_pct_val];

[cty_out_dist_val, cty_out_dist_ts, ~] = xlsread(cty_outage_dist_path);
cty_out_dist = [zeros(11*24, bus_county_num); cty_out_dist_val];

% fill in missing data
for col = 1:192
    for row = 2:455
        if cty_out_dist(row,col) == 0 && cty_out_dist(row-1, col) > 0 && sum(cty_out_dist(454:456, col)) > 0
            cty_out_dist(row,col) = cty_out_dist(row-1, col);
        end
    end
end

bus_out_dist = zeros([size(cty_out_dist,1), 2000]);

for c = 1:bus_county_num
    curr_cty = cty_out_header(c);
    bus_curr_cty = strcmp(bus_county, curr_cty);
    bus_cty_ratio = mpc.bus(bus_curr_cty, PD) / sum(mpc.bus(bus_curr_cty, PD));
    
    bus_out_dist(:, bus_curr_cty) = cty_out_dist(:, c) * bus_cty_ratio';
    
end
bus_out_dist(isnan(bus_out_dist)) = 0;

% computea area distribution using bus distribution
area_out_dist = zeros([size(cty_out_dist,1), 8]);
for a = 1:2000
    ar = mod(mpc.bus(a, BUS_AREA), 300);
    area_out_dist(:, ar) = area_out_dist(:, ar) + bus_out_dist(:, a);
end
area_out_dist = area_out_dist ./ sum(area_out_dist,2);
area_out_dist(isnan(area_out_dist)) = 0;

% outage data by county

[out_ng_val, out_ng_ts, ~] = xlsread(out_ng_path);
county_list = out_ng_ts(1,2:103);
county_num = size(county_list,2);
out_ng_val = [zeros(11*24, 2*county_num); out_ng_val];

[out_coal_val, out_coal_ts, ~] = xlsread(out_coal_path);
out_coal_val = [zeros(11*24, 2*county_num); out_coal_val];
[out_pv_val, out_pv_ts, ~] = xlsread(out_pv_path);
out_pv_val = [zeros(11*24, 2*county_num); out_pv_val];

[ng_derate_val, ng_derate_ts, ~] = xlsread(ng_derate_path);
ng_derate_val = [ones(11*24, county_num); ng_derate_val];

[~, gen_county, ~] = xlsread(gen_county_path);

dcnet_ts = lf_val(:,4);

wind_ts = readmatrix(wind_path);
wind_ts(:, 1) = [];
wind_ts(1, :) = [];
wind_ts(1:6,:) = []; % remove first 6 rows

solar_ts = readmatrix(solar_path);
solar_ts(:, 1) = [];
solar_ts(1, :) = [];

hydro_ts = readmatrix(hydro_path);
hydro_ts(:, 1) = [];
hydro_ts(1, :) = [];

genmix_header = ["wind", "solar", "hydro", "other", "ng", "coal", "nuclear"];

coal_ind = find(strcmp(mpc.genfuel, 'coal') == 1);
pv_ind = find(strcmp(mpc.genfuel, 'solar') == 1);
hydro_ind = find(strcmp(mpc.genfuel, 'hydro') == 1);
ng_ind = find(strcmp(mpc.genfuel, 'ng') == 1);
nuke_ind = find(strcmp(mpc.genfuel, 'nuclear') == 1);
wind_ind = find(strcmp(mpc.genfuel, 'wind') == 1);


load_rearranged = [load_val(:,3),load_val(:,4),load_val(:,8),load_val(:,6),load_val(:,5),load_val(:,7),load_val(:,1),load_val(:,2)];

% computer ercot load shed
ect_shed_ttl = lcf_ts(1:409,:) - load_rearranged(1:409,:);
ect_shed = ect_shed_ttl ./ sum(ect_shed_ttl,2);
ect_shed(isnan(ect_shed)) = 0;

% append non-load buses to population
all_load_idx = find(mpc.bus(:, PD) > 0);

bus_pop_full = zeros(bus_num, 4);
bus_pop_full(all_load_idx, :) = bus_pop(:, 2:5); 

bus_pop_area_pct = zeros(bus_num, 4);

% refactor population matrix
for load_area_i = 1:8
    real_area = load_area_i + 300;
    load_ind = find(mpc.bus(:, BUS_AREA) == real_area);
       
    
    for sec_i = 1:4
        pop_tmp = sum(bus_pop_full(load_ind, sec_i));
        bus_pop_area_pct(load_ind, sec_i) = bus_pop_full(load_ind, sec_i) / pop_tmp;
    end
end

% ercot load resources
bus_mw_total = zeros(504,1);
for load_area_i = 1:8
    sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
    bus_mw_total = bus_mw_total + sec_pct_all(1:504, 2) .* lcf_ts(:, load_area_i);
end
bus_mw_total(265:409) = bus_mw_total(265:409) + ers_ts(1:145) + rrs_ts(1:145);
clc;
case_path = "..\case\case_Mar3_5pm.mat";
mpc = loadcase(case_path);
 
%% modify base case gen
mpc.gen(:, GEN_STATUS) = 1;
gen_num = size(mpc.gen, 1);
% change the generator cost to linear (removing c2 in mpc.gencost)
mpc.gen(555, VG) = 1.02;

% save basic thermal capacity
gen_base = mpc.gen(:, PMAX);
gen_min = mpc.gen(:, PMIN);
coal_base = mpc.gen(coal_ind,PMAX);
coal_base(coal_base==0) = 1;
coal_min = mpc.gen(coal_ind,PMIN);
ng_base = mpc.gen(ng_ind,PMAX);
ng_base(ng_base==0) = 1;
ng_min = mpc.gen(ng_ind,PMIN);
nuke_base = mpc.gen(379, PMAX);
nuke_min = mpc.gen(379, PMIN);
wind_base = mpc.gen(wind_ind,MBASE);

% compute the area of generators
gen_area = floor(mod(mpc.gen(:, 1), 3000000) / 1000);

%% modify gen out data
out_ts(isnan(out_ts)) = 0;
%out_ts(:, 6) = out_ts(:, 6) + 500;
%out_ts(:, 5) = out_ts(:, 5) + 750

%out_ts(339, :) = out_ts(338, :);
%out_ts(340:342, :) = out_ts(341:343,:);
%out_ts(340:341, :) = out_ts(340:341, :) + 250;

%% modify DC tie flow
dcbus_ind = [145; 1966];
dc_ratio = [0.25, 0.75];



% merge line upgrades
%d3mpc = load('..\..\blackout\bte2k\hvdc_upgraded.mat');
%d3mpc = d3mpc.mpc_orig;
%mpc.branch(:,RATE_A) = max(mpc.branch(:,RATE_A), d3mpc.branch(d3mpc.branch(:, 1) > 3000000,RATE_A));
%mpc.branch(:,BR_X) = min(mpc.branch(:,BR_X), d3mpc.branch(d3mpc.branch(:, 1) > 3000000,BR_X));


%% DR related variables and containers

% load rationing
%lr_pct_max = 0.1;
lr_pct_step = 0.005;
%lr_pct = 0;


% ERCOT load resources
ercot_lr_enable = 1;
%ercot_lr_scaler = 1;

ercot_lr_pct_max = 1;
ercot_lr_pct_step = 0.05;
ercot_lr_max_ramp_up = ercot_lr_pct_max;
ercot_lr_max_ramp_dn = ercot_lr_pct_max / 2;

% incentive based
inc_enable = 1;
%res_participation_ratio = 0.1; % maximum proportion of res participants

%inc_bias = 0.85; % fix the randomization to be the m%ean +- std;%
%inc_bias = 1.15;
%inc_bias = 1;

active_ratio = 0.25;


inc_rnd_std_pct = 0; % percentage of mean



% loop to create more data
for  lr_pct_max= 0.50:-0.05:0.4
    for res_participation_ratio = 0.5:-0.05:0.15
        for ercot_lr_scaler = 6.5:-0.3:1
 %           for inc_bias = [1]
                inc_bias = 1;
                clc;

                lr_pct = 0;
                lr_max_ramp_up = lr_pct_max / 6;
                lr_max_ramp_dn = lr_pct_max / 6;
                

                ercot_lr_pct = 0;
                lr_shed_total = zeros(bus_num, 1);
                inc_shed_act = zeros(bus_num, 1);
                inc_shed_inact = zeros(bus_num, 1);
                
                act_red_pcts = inc_bias * [0.726639835000000;0.702821616000000;0.717001661000000;0.836879948000000;0.763512598000000;0.751789263000000;0.884280199000000;0.906001582000000;0.912319039000000;0.939109041000000;0.940497922000000;0.929441939000000;0.931010453000000;0.813742696000000;0.796455138000000;0.777696394000000;0.681431371000000;0.678340452000000;0.741335527000000;0.665097755000000;0.647544457000000;0.638205616000000;0.655924266000000;0.657043196000000];
                inact_red_pcts = inc_bias * [0.940589702000000;0.947639628000000;0.939089080000000;0.913582717000000;0.939623458000000;1;1;0.944609072000000;1;1;1;1;1;1;1;0.964302885000000;1;1;0.939092357000000;0.942080133000000;0.915094340000000;0.829204954000000;0.852977879000000;0.859254734000000];

                
                ercot_available_mw = min((ers_ts(1:145) + rrs_ts(1:145))*ercot_lr_scaler, bus_mw_total(265:409));
                
                outfolder = sprintf('../border/%2.2f_%2.2f_%2.2f', lr_pct_max, ercot_lr_scaler*ercot_lr_enable, res_participation_ratio*inc_enable);
                
                disp(outfolder);
                %% looping through hours
                eta = 0;
                % Start with 2/1 0:00
                total_rows = size(load_val, 1);
                real_date = 1;
                real_hr = 0;

                % Start with 2/12 0:00
                day_tar = 12; % Feb 12
                hour_tar = 0;
                row_num_start = (day_tar-1) * 24 + hour_tar + 1;
                time_steps = (row_num_start:409)'; %load data only up to 2/18

                % initialize containers
                ts_num = size(time_steps,1);
                system_gen = zeros(ts_num,1);
                system_load = zeros(ts_num,1);
                genmix_mismatch = zeros(ts_num, 7);
                genmix_opf = zeros(ts_num, 7);
                genmix_cf = zeros(ts_num, 7);
                loaddiff_cf = zeros(ts_num, 1);

                all_gen_pg = zeros(ts_num, size(mpc.gen, 1));
                all_brn_pf = zeros(ts_num, size(mpc.branch, 1));


                % DR containers
                lr_pcts = zeros(ts_num, 1);
                ercot_lr_pcts = zeros(ts_num, 1);
                ci = 1; % container index

                lr_mws = zeros(ts_num, 1);
                ercot_lr_mws = zeros(ts_num, 1);

                total_shed_mws = zeros(ts_num, 1);
                inc_mws = zeros(ts_num, 1);

                % load shed container
                real_shed = zeros(ts_num, 8);
                forced_shed = zeros(ts_num, 1);
                area_ratios = zeros(ts_num, 8);

                for row_num = row_num_start:409
                    % fetch and match gen
                    % reset gen PMIN
                    %disp(row_num);

                    mpc.gen(1:gen_num, PMIN) = gen_min;

                    for gen_type_i = 1:7
                        gen_type = genmix_header(gen_type_i);
                        if strcmp(gen_type,'other')

                        else
                            type_ind = find(strcmp(mpc.genfuel, gen_type) == 1);

                            % set pmin/max if renewable(allow curtailment), just set PMAX
                            if gen_type_i == 1 % wind
                                % scale wind according to real EIA data AND anti-frost
                                bes_wind = sum(wind_ts(row_num,:));
                                eia_wind = gen_val(row_num, gen_type_i);
                                real_wind_pct = eia_wind / bes_wind;
                                mpc.gen(type_ind, PG) = wind_ts(row_num,:)' .* real_wind_pct;
                                mpc.gen(type_ind, PMIN) = 0;
                                mpc.gen(type_ind, PMAX) = mpc.gen(type_ind, PG);

                            elseif gen_type_i == 2 % solar
                                % scale PV according to real EIA
                                bes_pv = sum(solar_ts(row_num,:));
                                eia_pv = gen_val(row_num, gen_type_i);
                                if bes_pv == 0
                                    bes_pv = 1;
                                end
                                real_pv_pct = eia_pv / bes_pv;
                                cf_pv_pct = 1 - (1 - real_pv_pct);
                                mpc.gen(type_ind, PG) = solar_ts(row_num,:)';
                                mpc.gen(type_ind, PMIN) = 0;
                                mpc.gen(type_ind, PMAX) = mpc.gen(type_ind, PG);   

                            elseif gen_type_i == 3 % hydro
                                mpc.gen(type_ind, PG) = hydro_ts(row_num,:)';
                                mpc.gen(type_ind, PMIN) = 0;
                                mpc.gen(type_ind, PMAX) = mpc.gen(type_ind, PG);

                            elseif gen_type_i == 5 % ng
                                for cty = 1:county_num
                                    % check capacity in cty
                                    ng_ttl = out_ng_val(row_num, cty+county_num);
                                    if ng_ttl == 0
                                        continue
                                    end

                                    % find all gens within the county
                                    curr_cty = county_list(cty);
                                    curr_gen_ind = intersect(find(strcmp(gen_county, curr_cty) == 1), type_ind);

                                    ng_out = out_ng_val(row_num, cty) * ng_derate_val(row_num,cty);

                                    % scale corresponding gens
                                    fac1 = (1 - ng_out / ng_ttl);
                                    fac2 = (sum(gen_base(curr_gen_ind)) - ng_out) / sum(gen_base(curr_gen_ind));
                                    fac = max([fac1, fac2]);
                                    mpc.gen(curr_gen_ind, PMAX) = gen_base(curr_gen_ind) .* fac;
                                    mpc.gen(curr_gen_ind, PMIN) = gen_min(curr_gen_ind) .* fac;
                                    % record total outed capacity in the specified area

                                end
                                % scale remaining gens using outage nums
                                ng_up = (sum(ng_base) - out_ts(row_num, gen_type_i));
                                mpc.gen(type_ind, PMAX) = mpc.gen(type_ind, PMAX) * ng_up / sum(mpc.gen(type_ind, PMAX));
                                % scale down PMIN as result of outage
                                mpc.gen(type_ind, PMIN) = mpc.gen(type_ind, PMAX) .* ng_min ./ ng_base;

                            elseif gen_type_i == 6 % coal
                                % gen outage by county
                                for cty = 1:county_num
                                    % check capacity in cty
                                    coal_ttl = out_coal_val(row_num, cty+county_num);
                                    if coal_ttl == 0
                                        continue
                                    end
                                    % find all gens within the county
                                    curr_cty = county_list(cty);
                                    curr_gen_ind = intersect(find(strcmp(gen_county, curr_cty) == 1), type_ind);

                                    coal_out = out_coal_val(row_num, cty);
                                    % scale corresponding gens
                                    fac1 = (1 - coal_out / coal_ttl);
                                    fac2 = (sum(gen_base(curr_gen_ind)) - coal_out) / sum(gen_base(curr_gen_ind));
                                    fac = max([fac1, fac2]);
                                    mpc.gen(curr_gen_ind, PMAX) = gen_base(curr_gen_ind) * fac;
                                    mpc.gen(curr_gen_ind, PMIN) = gen_min(curr_gen_ind) * fac;
                                end
                                % scale remaining gens using outage nums
                                coal_up = (sum(coal_base) - out_ts(row_num, gen_type_i));
                                mpc.gen(type_ind, PMAX) = mpc.gen(type_ind, PMAX) * coal_up / sum(mpc.gen(type_ind, PMAX));
                                mpc.gen(type_ind, PMIN) = mpc.gen(type_ind, PMAX) .* coal_min ./ coal_base;
                                % fix PMAX < PMIN

                            elseif gen_type_i == 7 % nuke
                                if out_ts(row_num, gen_type_i) > nuke_base
                                    mpc.gen(379, GEN_STATUS) = 0;
                                    mpc.gen(379, PMAX) = 0;
                                else
                                    mpc.gen(379, GEN_STATUS) = 1;
                                    mpc.gen(379, PMAX) = nuke_base - out_ts(row_num, gen_type_i);
                                    mpc.gen(379, PMIN) = nuke_min - out_ts(row_num, gen_type_i);
                                end
                            end
                        end
                    end

                    % reset existing DC injection
                    %mpc.bus(dcbus_ind, PD) = 0;

                    % fetch and match load
                    for load_area_i = 1:8
                        real_area = load_area_i + 300;
                        load_ind = find(mpc.bus(:, BUS_AREA) == real_area);
                        base_load = sum(mpc.bus(load_ind, PD));

                        % load shed

                        % scale load shed inherited from last step

                        if ci > 1
                            real_shed(ci, load_area_i) = real_shed(ci - 1, load_area_i) * lcf_ts(row_num, load_area_i) / lcf_ts(row_num-1, load_area_i);
                        end
                        real_load = lcf_ts(row_num, load_area_i) - real_shed(ci, load_area_i);

                        % scale load
                        mpc.bus(load_ind, PD) = mpc.bus(load_ind, PD) * real_load / base_load;
                        %mpc.bus(load_ind, PD) = all_bus_pd(load_ind) * real_load / base_load;

                        % distribute load among buses based on population
                        %area_pop = sum(bus_pop_full(load_ind, 4));
                        %mpc.bus(load_ind, PD) = bus_pop_full(load_ind, 4) * real_load / area_pop;
                    end

                    % (TEST) 123333333333333333333333333333333333333
                    %res_tot_mw = 0;
                    %for load_area_i = 1:8
                    %    sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                    %    res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);
                    %    res_tot_mw = res_tot_mw + res_mw;
                    %end
                    %
                    %lr_pct_max = ercot_available_mw(ci) / res_tot_mw;
                    %lr_pct = min(lr_pct_max, lr_pct);
                    
                    
                    % load including forced shedding from the past
                    load_before_lr = mpc.bus(:, PD);

                    % compute inherited load reduction and reduce from LR
                    res_tot_mw = 0;
                    for load_area_i = 1:8
                        real_area = load_area_i + 300;
                        load_ind = find(mpc.bus(:, BUS_AREA) == real_area);

                        % load rationing
                        sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                        res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);
                        res_tot_mw = res_tot_mw + res_mw;

                        % pop dist
                        %lr_shed_total(load_ind) = lr_pct * res_mw * bus_pop_area_pct(load_ind, 1);

                        % pd dist
                        lr_shed_total(load_ind) = lr_pct * res_mw * load_before_lr(load_ind) / sum(load_before_lr(load_ind));
                    end
                    
                    mpc.bus(:, PD) = load_before_lr - lr_shed_total;


                    % inherited ERCOT load resources
                    if ercot_lr_enable

                        % total business MW from all areas

                        ercot_lr_mw = ercot_available_mw(ci) * ercot_lr_pct;


                        for load_area_i = 1:8
                            load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));

                            sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                            bus_mw_pct = sec_pct_all(row_num, 2) * lcf_ts(row_num, load_area_i) / bus_mw_total(row_num);

                            ercot_lr_area = ercot_lr_mw * bus_mw_pct;
                            mpc.bus(load_ind, PD) = mpc.bus(load_ind, PD) * (1 - ercot_lr_area / sum(mpc.bus(load_ind, PD)));
                        end

                    end

                    % modify existing DC injection
                    %mpc.bus(dcbus_ind, PD) = dc_ratio * dcnet_ts(row_num);

                    % add DC as generators
                    mpc.gen = [mpc.gen; zeros(2,25)];
                    mpc.gen(607,1) = mpc.bus(dcbus_ind(1), 1);
                    mpc.gen(608,1) = mpc.bus(dcbus_ind(2), 1);
                    mpc.gen([607, 608],PMAX) = dc_ratio * -dcnet_ts(row_num);
                    mpc.gen([607, 608],PMIN) = dc_ratio * -dcnet_ts(row_num);
                    mpc.gen(607,19) = inf;
                    mpc.gen(608,19) = inf;
                    mpc.gen(607,GEN_STATUS) = 1;
                    mpc.gen(608,GEN_STATUS) = 1;


                    mpc.genfuel = [mpc.genfuel; "hvdc"; "hvdc"];
                    mpc.gencost = [mpc.gencost; [2,0,0,3,0,0,0]; [2,0,0,3,0,0,0]];   

                    % solve OPF 
                    opfres = rundcopf(mpc, mpopt);
                    %assert(opfres.success);

                    if row_num == 398
                        disp(123);
                    end

                    % if opf fails, allocate DR
                    ercot_lr_pct_old = ercot_lr_pct;
                    % ERCOT load resources
                    if ~opfres.success && ercot_lr_pct < ercot_lr_pct_max && ercot_lr_enable

                        while ~opfres.success && ercot_lr_pct < ercot_lr_pct_max && ercot_lr_pct - ercot_lr_pct_old < ercot_lr_max_ramp_up


                            ercot_lr_pct_diff = min(ercot_lr_pct + ercot_lr_pct_step, ercot_lr_pct_max) - ercot_lr_pct;
                            ercot_lr_pct = min(ercot_lr_pct + ercot_lr_pct_step, ercot_lr_pct_max);

                            ercot_lr_mw_diff = ercot_available_mw(ci) * ercot_lr_pct_diff;

                            %mpc.bus(:, PD) = mpc.bus(:, PD) * (1 - ercot_lr_mw_diff / sum(mpc.bus(:, PD)));

                            for load_area_i = 1:8
                                load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));

                                sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                                bus_mw_pct = sec_pct_all(row_num, 2) * lcf_ts(row_num, load_area_i) / bus_mw_total(row_num);

                                ercot_lr_area_diff = ercot_lr_mw_diff * bus_mw_pct;
                                mpc.bus(load_ind, PD) = mpc.bus(load_ind, PD) * (1 - ercot_lr_area_diff / sum(mpc.bus(load_ind, PD)));
                            end

                            opfres = rundcopf(mpc, mpopt);
                        end
                    end
                    ercot_lr_pcts(ci) = ercot_lr_pct;

                    % incentive based
                    hour = mod((ci-1), 24) + 1;
                    inc_shed_total = zeros(bus_num, 1);
                    all_trend = 1 - 0.2 * sum(real_shed(ci, :)) / 25000;
                    curr_reserve = sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD));

                    if (~opfres.success || sum(real_shed(ci, :)) > 10 || curr_reserve<2000) && inc_enable
                        curr_reserve = max(sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD)), 0);
                        % get corresponding residential load to respond randomly
                        for load_area_i = 1:8
                            load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));

                            % available MW for inc
                            sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                            res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);
                            inc_tar_mw = min(res_mw * res_participation_ratio * all_trend, max(2000-curr_reserve*lcf_ts(row_num, load_area_i)/sum(lcf_ts(ci,:)),0));

                            act_real_mw = inc_tar_mw * active_ratio;
                            inact_real_mw = inc_tar_mw * (1 - active_ratio);

                            % pd based
                            inc_shed_act(load_ind) = act_real_mw * mpc.bus(load_ind, PD) / sum(mpc.bus(load_ind, PD));
                            inc_shed_inact(load_ind) = inact_real_mw * mpc.bus(load_ind, PD) / sum(mpc.bus(load_ind, PD));

                            % pop based
                            %inc_shed_act(load_ind) = act_real_mw * bus_pop_area_pct(load_ind, 1);
                            %inc_shed_inact(load_ind) = inact_real_mw * bus_pop_area_pct(load_ind, 1);
                        end

                        % random process
                        inc_shed_act = inc_shed_act .* normrnd(act_red_pcts(hour), act_red_pcts(hour)*inc_rnd_std_pct, bus_num, 1);
                        inc_shed_inact = inc_shed_inact .* normrnd(inact_red_pcts(hour), inact_red_pcts(hour)*inc_rnd_std_pct, bus_num, 1);

                        inc_shed_total = inc_shed_act + inc_shed_inact;

                        % check if lr > bus PD
                        inc_shed_total = min(mpc.bus(:, PD), inc_shed_total);

                        % apply to pd
                        mpc.bus(:, PD) = mpc.bus(:, PD) - inc_shed_total;

                    end
                    inc_mws(ci) = sum(inc_shed_total);
                    opfres = rundcopf(mpc, mpopt);

                    % load including forced shedding from the past
                    load_before_lr = mpc.bus(:, PD);

                    % load rationing
                    lr_pct_old = lr_pct;
                    lr_shed_total = zeros(bus_num,1);
                    if ~opfres.success && lr_pct < lr_pct_max
                        while ~opfres.success && lr_pct < lr_pct_max && (lr_pct - lr_pct_old) < lr_max_ramp_up
                            lr_pct = min(lr_pct + lr_pct_step, lr_pct_max);

                            % simple iteration through areas
                            res_tot_mw = 0;
                            for load_area_i = 1:8
                                load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));

                                % load rationing
                                sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                                res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);
                                res_tot_mw = res_tot_mw + res_mw;

                                %lr_shed_total(load_ind) = lr_pct * res_mw *  bus_pop_area_pct(load_ind, 1);
                                lr_shed_total(load_ind) = (lr_pct - lr_pct_old) * res_mw * load_before_lr(load_ind) / sum(load_before_lr(load_ind));
                            end
                            % check if lr > bus PD
                            %lr_shed_total = min(load_before_lr, lr_shed_total);
                            %lr_shed_total(dcbus_ind) = 0;

                            % allocate remaining to all buses
                            %lr_residue = lr_pct * res_tot_mw - sum(lr_shed_total);
                            %rem_loads = load_before_lr - lr_shed_total;
                            %lr_shed_total = lr_shed_total + lr_residue * rem_loads / sum(rem_loads);

                            % apply!
                            mpc.bus(:, PD) = load_before_lr - lr_shed_total;

                            % try opf
                            opfres = rundcopf(mpc, mpopt);
                        end
                    end



                    % calculate area base load for load shed
                    area_ratio = zeros(1,8);
                    area_load_base = zeros(1,8);
                    for load_area_i = 1:8
                        load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));
                        area_load_base(load_area_i) = sum(mpc.bus(load_ind, PD));
                        area_ratio(load_area_i) = sum(mpc.bus(load_ind, PD)) / sum(mpc.bus(:, PD));
                    end

                    area_ratio = (1 - eta)*area_ratio + eta*ect_shed(row_num, :);
                    area_ratios(ci, :) = area_ratio;


                    % if opf fails after DR, shed load
                    if ~opfres.success
                       %disp('error'); 
                       %disp(row_num);
                       % % use storage
                       % mpc.bus(mpc.bus(:,PD) > 0, PD) = mpc.bus(mpc.bus(:,PD) > 0, PD) * (sum(mpc.bus(mpc.bus(:,PD) > 0, PD)) - 1000) / sum(mpc.bus(mpc.bus(:,PD) > 0, PD));

                       % add load shed (by area)
                       % start from capacity needed
                        base_shed = max(sum(mpc.bus(:,PD)) - sum(mpc.gen(:,PMAX)), 0); 
                        load_before_shed = mpc.bus(:, PD);
                        shed_step = 100; % MW
                        new_shed = zeros(1,8);

                        % assign base shed load to the ratio
                        new_shed(1, :) = area_ratio .* base_shed;
                        while ~opfres.success 
                            new_shed(1, :) = new_shed(1, :) + area_ratio .* shed_step;
                            % simple iteration through areas
                            for load_area_i = 1:8
                                load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));
                                shed_fac = (area_load_base(load_area_i) - new_shed(1, load_area_i)) / area_load_base(load_area_i);
                                mpc.bus(load_ind, PD) = load_before_shed(load_ind) * shed_fac;
                            end
                            % try opf
                            opfres = rundcopf(mpc, mpopt);
                        end
                        real_shed(ci, :) = real_shed(ci, :) + new_shed;
                    % if success, check if can reconnect
                    elseif row_num > 340
                        reconn_step = 100; % MW
                        min_reserve = 000; % minimum reserve
                        max_reconn = 1000; % maximum reconnection per hr

                        curr_reserve = max(sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD)), 0);
                        load_before_reconn = mpc.bus(:, PD);
                        % start reconnecting if margin > 1000
                        if curr_reserve - min_reserve > 3000 && sum(real_shed(ci,:)) > 0
                            % reconnect iteratively 
                            new_reconn = zeros(1,8);
                            while opfres.success && sum(new_reconn) < 1000 && sum(new_reconn) < sum(real_shed(ci,:))
                                new_reconn(1, :) = new_reconn(1, :) + area_ratio .* reconn_step;
                                if sum(new_reconn(1,:)) > sum(real_shed(ci,:))
                                    new_reconn(1,:) = real_shed(ci,:);
                                end
                                for load_area_i = 1:8
                                    if new_reconn(1,load_area_i) > real_shed(ci, load_area_i)
                                        new_reconn(1,load_area_i) = real_shed(ci, load_area_i);
                                    end 
                                end
                                % simple iteration through areas
                                for load_area_i = 1:8
                                    load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));
                                    reconn_fac = (area_load_base(load_area_i) + new_reconn(1, load_area_i)) / area_load_base(load_area_i);
                                    mpc.bus(load_ind, PD) = load_before_reconn(load_ind) * reconn_fac;
                                end
                                % try opf
                                opfres = rundcopf(mpc, mpopt);               
                            end
                            real_shed(ci, :) = real_shed(ci, :) - new_reconn;
                        end  
                    end 


                    % try reduce load rationing when there is no forced load shed
                %     if opfres.success && lr_pct > 0 && (sum(real_shed(ci, :)) < 10)
                %         
                %         while opfres.success && lr_pct > 0 
                %             lr_pct = max(lr_pct - lr_pct_step, 0);
                %             
                %             % simple iteration through areas
                %             for load_area_i = 1:8
                %                 load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));
                % 
                %                 % load rationing
                %                 sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                %                 res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);
                %                 
                %                 %lr_shed_total(load_ind) = lr_pct * res_mw *  bus_pop_area_pct(load_ind, 1);
                %                 lr_shed_total(load_ind) = lr_pct * res_mw * load_before_lr(load_ind) / sum(load_before_lr(load_ind));
                %             end
                %             mpc.bus(:, PD) = load_before_lr - lr_shed_total;
                %             
                %             % try opf
                %             opfres = rundcopf(mpc, mpopt);            
                %         end
                %     end
                    curr_reserve = max(sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD)), 0);

                    lr_pct_old = lr_pct;
                    if opfres.success && lr_pct > 0 && (sum(real_shed(ci, :)) < 10) && curr_reserve > 2000

                        while opfres.success && lr_pct > 0 && (lr_pct_old - lr_pct) < lr_max_ramp_dn && curr_reserve > 2000
                            % save pre-reconn vars
                            pd_tmp = mpc.bus(:, PD);
                            lr_shed_total_tmp = lr_shed_total;
                            lr_pct_tmp = lr_pct;

                            lr_pct_diff = lr_pct - max(lr_pct - lr_pct_step, 0);

                            % simple iteration through areas
                            for load_area_i = 1:8
                                load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));

                                % load rationing
                                sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                                res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);

                                %lr_shed_total(load_ind) = lr_pct * res_mw * load_before_lr(load_ind) / sum(load_before_lr(load_ind));
                                %lr_shed_total(load_ind) = lr_shed_total(load_ind) * (lr_pct - lr_pct_diff) / lr_pct;
                                lr_reconn_mw = res_mw * lr_pct_diff * mpc.bus(load_ind, PD) / sum(mpc.bus(load_ind, PD));
                                lr_shed_total(load_ind) = lr_shed_total(load_ind) + lr_reconn_mw;
                                mpc.bus(load_ind, PD) = mpc.bus(load_ind, PD) + lr_reconn_mw;
                            end

                            lr_pct = max(lr_pct - lr_pct_step, 0);

                            % try opf
                            opfres = rundcopf(mpc, mpopt);            
                            curr_reserve = max(sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD)), 0);

                        end

                        % reverse 1 step if opf fails
                        if ~opfres.success
                            lr_pct = lr_pct_tmp;
                            lr_shed_total = lr_shed_total_tmp;
                            mpc.bus(:, PD) = pd_tmp;

                            opfres = rundcopf(mpc, mpopt);  
                            %assert(opfres.success);
                        end
                    end

                    % reduce ERCOT load resource shed
                    curr_reserve = max(sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD)), 0);

                    ercot_lr_pct_old = ercot_lr_pct;
                    if opfres.success && ercot_lr_pct > 0 && lr_pct < 0.01 && ercot_lr_enable && curr_reserve > 3000 && (sum(real_shed(ci, :)) < 10)

                        while opfres.success && ercot_lr_pct > 0 && (ercot_lr_pct_old - ercot_lr_pct) < ercot_lr_max_ramp_dn && curr_reserve > 3000
                            % store pre-reconn vars
                            pd_tmp = mpc.bus(:, PD);
                            ercot_lr_pct_tmp = ercot_lr_pct;

                            ercot_lr_pct_diff = ercot_lr_pct - max(ercot_lr_pct - ercot_lr_pct_step, 0);

                            ercot_lr_mw_reconn = ercot_available_mw(ci) * ercot_lr_pct_diff;
                            %mpc.bus(:, PD) = mpc.bus(:, PD) * (1 + ercot_lr_mw_reconn / sum(mpc.bus(:, PD)));


                            for load_area_i = 1:8
                                load_ind = find(mpc.bus(:, BUS_AREA) == (load_area_i + 300));

                                sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                                bus_mw_pct = sec_pct_all(row_num, 2) * lcf_ts(row_num, load_area_i) / bus_mw_total(row_num);

                                ercot_lr_area_reconn = ercot_lr_mw_reconn* bus_mw_pct;
                                mpc.bus(load_ind, PD) = mpc.bus(load_ind, PD) * (1 + ercot_lr_area_reconn / sum(mpc.bus(load_ind, PD)));
                            end


                            ercot_lr_pct = max(ercot_lr_pct - ercot_lr_pct_step, 0);
                            opfres = rundcopf(mpc, mpopt);
                            curr_reserve = max(sum(mpc.gen(:,PMAX)) - sum(mpc.bus(:,PD)), 0);
                        end

                        % reverse 1 step if opf fails
                        if ~opfres.success
                            mpc.bus(:, PD) = pd_tmp;
                            ercot_lr_pct = ercot_lr_pct_tmp;
                            opfres = rundcopf(mpc, mpopt);  
                            assert(opfres.success);
                        end
                    end


                    lr_pcts(ci) = lr_pct;
                    ercot_lr_pcts(ci) = ercot_lr_pct;

                    % remove dc gens
                    mpc.gen([607;608],:) = [];
                    mpc.genfuel([607;608],:) = [];
                    mpc.gencost([607;608],:) = [];




                    % SAVE RESULTS
                    % calculate each type of addition
                    %if row_num > 313    
                    all_gen_pg(ci, :) = opfres.gen(1:gen_num, PG)';
                    all_brn_pf(ci, :) = opfres.branch(:, PF)';

                    % plot load and demand
                    system_gen(ci, 1) = sum(opfres.gen([coal_ind;ng_ind;nuke_ind],PMAX)) + sum(opfres.gen([wind_ind;pv_ind;hydro_ind],PG));
                    system_load(ci, 1) = sum(opfres.bus(:,PD));

                    % demand response related
                    % load rationing
                    res_tot_mw = 0;
                    for load_area_i = 1:8
                        sec_pct_all = area_profiles_pct(area_mpc(load_area_i));
                        res_mw = sec_pct_all(row_num, 1) * lcf_ts(row_num, load_area_i);
                        res_tot_mw = res_tot_mw + res_mw;
                    end
                    lr_mws(ci) = res_tot_mw * lr_pct;

                    % ercot load resources
                    ercot_lr_mws(ci) = ercot_lr_pct * ercot_available_mw(ci);

                    total_shed_mws(ci) = sum(real_shed(ci, :)) + lr_mws(ci) + ercot_lr_mws(ci) + inc_mws(ci);
                    forced_shed(ci, 1) = sum(real_shed(ci, :));


                    %total_shed_mws(ci) = sum(lcf_ts(row_num, :)) - system_load(ci, 1);
                    %forced_shed(ci, 1) = total_shed_mws(ci) - lr_mws(ci) - ercot_lr_mws(ci) - inc_mws(ci);

                    % gen difference by type
                    for gen_type_i = 1:7
                        gen_type = genmix_header(gen_type_i);
                        if gen_type == 'other'
                            ;
                        else
                            type_ind = find(strcmp(mpc.genfuel, gen_type) == 1);
                            real_gen = gen_val(row_num, gen_type_i);
                            opf_gen = sum(opfres.gen(type_ind, PG));
                            genmix_mismatch(ci, gen_type_i) = opf_gen - real_gen;
                            genmix_opf(ci, gen_type_i) = opf_gen;
                            % save counterfactural
                            genmix_cf(ci, gen_type_i) = opf_gen - real_gen;
                        end
                    end
                    loaddiff_cf(ci, 1) = sum(lcf_ts(row_num, :)) - sum(load_rearranged(row_num, :));

                    real_hr = real_hr + 1;
                    ci = ci + 1;
                    if real_hr == 24
                        real_hr = 0;
                        real_date = real_date + 1;
                    end


                end

                %all_gen_pg = [double(mpc.genid'); all_gen_pg];
                %all_brn_pf = [double(mpc.branchid'); all_brn_pf];
                system_gen(1:48) = system_gen(1:48) + (system_gen(49) - system_gen(48)); % fix offset before Feb. 14 (lack of real data)
                ercot_load = sum(load_val,2);

                final_shed = sum(real_shed, 2);

                %% output csvs
                
                if sum(forced_shed > 100)
                
                    mkdir(outfolder);

                    stamps = gen_ts(time_steps+24, :);
                    total_shed_out = [stamps, num2cell(total_shed_mws)];
                    forced_shed_out = [stamps, num2cell(forced_shed)];
                    load_rationing_pcts_out =[stamps, num2cell(lr_pcts)];
                    load_rationing_mws_out = [stamps, num2cell(lr_mws)]; 
                    ercot_lr_pcts_out = [stamps, num2cell(ercot_lr_pcts)];
                    ercot_lr_mws_out = [stamps, num2cell(ercot_lr_mws)];
                    inc_mws_out = [stamps, num2cell(inc_mws)];

                    writecell(total_shed_out, outfolder + "/total_shed_mws.csv");
                    %writecell(forced_shed_out, outfolder + sprintf("/forced_shed_mws_%2.2f.csv", inc_bias));
                    writecell(forced_shed_out, outfolder + sprintf("/forced_shed_mws.csv"));
                    writecell(load_rationing_pcts_out, outfolder + "/load_rationing_pcts.csv");
                    writecell(load_rationing_mws_out, outfolder + "/load_rationing_mws.csv");
                    writecell(ercot_lr_pcts_out, outfolder + "/ercot_load_resources_pcts.csv");
                    writecell(ercot_lr_mws_out, outfolder + "/ercot_load_resources_mws.csv");
                    %writecell(inc_mws_out, outfolder + sprintf("/incentivized_mws_%2.2f.csv", inc_bias));
                    writecell(inc_mws_out, outfolder + sprintf("/incentivized_mws.csv"));
                    
                    break;
                end
%            end
        end
    end
end
%lcf_total = sum(lcf_ts, 2);
%lcf_total = lcf_total(time_steps);
%lcf_total - system_load - total_shed_mws;
%% plot results

%plot summary curves
% ax1 = subplot(3,1,1);
% plot(ax1,time_steps, system_gen);
% hold on;
% plot(ax1,time_steps, system_load);
% plot(ax1,time_steps, cap_ts(time_steps,1));
% plot(ax1,time_steps, ercot_load(time_steps));
% title(ax1, 'Simulated Capacity and Load During the Event');
% xlabel(ax1, 'Time / hr');
% ylabel(ax1, 'Capacity / MW');
% legend(ax1, ["Gen","Load","ERCOT-Capacity","ERCOT-Load"]);
% 
% %ax2 = subplot(3,1,2);
% %plot(ax2, time_steps, genmix_mismatch);
% %title(ax2, 'Dispatch Mismatch grouped by Generation Source');
% %xlabel(ax2, 'Time / hr');
% %ylabel(ax2, 'Mismatch / MW');
% %legend(ax2, genmix_header);
% 
% ax2 = subplot(3,1,2);
% plot(ax2, time_steps, total_shed_mws);
% hold on;
% plot(ax2, time_steps, lr_mws);
% plot(ax2, time_steps, ercot_lr_mws);
% plot(ax2, time_steps, inc_mws);
% plot(ax2, time_steps, forced_shed);
% title(ax2, 'Decomposition of Load Shedding');
% xlabel(ax2, 'Time / hr');
% ylabel(ax2, 'Capacity / MW');
% legend(ax2, ["Total","Load Rationing", "ERCOT Load Resources", "Incentive Based", "Forced Shed"]);
% 
% ax3 = subplot(3,1,3);
% plot(ax3, time_steps, total_shed_mws);
% hold on;
% plot(ax3, time_steps, shed_ts(time_steps, 2));
% title(ax3, 'Real and Simulated Load Shed');
% xlabel(ax3, 'Time / hr');
% ylabel(ax3, 'Capacity / MW');
% legend(ax3, ["Simulated","Real"]);
% 

% 
% % plot bar chart for load shedding
% lcf_ttl = sum(lcf_ts(time_steps,:) ,2);
% real_shed_nml = area_ratios(77:end,:);
% ect_shed_nml = ect_shed_ttl(time_steps,:) ./ sum(ect_shed_ttl(time_steps,:),2);
% ect_shed_nml = ect_shed_nml(77:end,:);
% 
% real_shed_means = mean(real_shed_nml, 1);
% real_shed_stds = 3 * std(real_shed_nml, 0, 1);
% ect_shed_means = mean(ect_shed_nml, 1);
% ect_shed_stds = 3 * std(ect_shed_nml, 0, 1);
% 
% figure();
% errorbar(1:8, real_shed_means, real_shed_stds,'-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on;
% errorbar(1:8, ect_shed_means, ect_shed_stds,'-s','MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
% 
% fig3b.ercot_shed = ect_shed_ttl(time_steps,:);
% fig3b.sim_shed = real_shed;
% fig3b.sim_load = lcf_ts(time_steps,:) - real_shed;
% fig3b.ercot_load = load_rearranged(time_steps,:);
% % 
% ercot_shed = fig3b.ercot_shed(72+20:72+36,:);
% ercot_load = fig3b.ercot_load(72+20:72+36,:);
% sim_shed = fig3b.sim_shed(72+20:72+36,:);
% sim_load = fig3b.sim_load(72+20:72+36,:);
% 
% ercot_shed_pct = ercot_shed./sum(ercot_shed,2);
% ercot_load_pct = ercot_load./sum(ercot_load,2);
% sim_shed_pct = sim_shed./sum(sim_shed,2);
% sim_load_pct = sim_load./sum(sim_load,2);
% 
% sim_odds_ratio = sim_shed_pct./sim_load_pct;
% ercot_odds_ratio = ercot_shed_pct./ercot_load_pct;
% 
% sim_odds_ratio_mean = mean(sim_odds_ratio);
% sim_odds_ratio_std = std(sim_odds_ratio);
% ercot_odds_ratio_mean = mean(ercot_odds_ratio); 
% ercot_odds_ratio_std = std(ercot_odds_ratio);
% x=1:8;
% 
% figure(4)
% bar(x, [ercot_odds_ratio_mean', sim_odds_ratio_mean']);