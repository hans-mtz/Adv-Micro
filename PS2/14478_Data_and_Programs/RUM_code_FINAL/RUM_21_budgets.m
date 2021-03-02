function [budgets,shares,income,Z] = RUM_21_budgets
%% Code Description: Create budges and expenditure shares
% This code loads the UK household-survey dataset.  It cleans the data by
% restricting to households with a car and at least one child.  
% Next, it computes expenditure shares for each of our classes (food, 
% non-durables, and services for the classes=3 case). 
% Last, it computes normalized budgets for each year using median
% expenditure.
% Output:
%   - budgets is a [periods, classes] matrix of normalized budgets
%   - shares is a [num_household,classes]{periods} structure of expenditure
%     shares 
%   - income is a [num_household,1]{periods} structure of household
%    expenditures.
%   - Z is the instrument used in endogenous series estimator.
%       It is total household income (not to be confused with the income
%       vector which is actually expenditure). 

%% Global parameters
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores

%% Classes and Catagories
% Below are column numbers for each catagory of expenditure.  These are
% used to aggregate into our classes.
% Example: 4 classes
    % -Food items: groceries 01-24, eat out 25-26
    % -clothing and shoes, 50-54
    % -other nondurables: alcohol and tobacco products 27-30, pets 44-45,
    %   chemist's good 55, leisure 64-67  
    % -services: fuels 36-39, communication 46-49, personal services 56, 
    % entertainment 68-69, insurance, maintainance 58-60, travel 61-63

% Code currently allows classes to be set to 3,4,5.  One can adjust the
% class and catagory definition below to get other aggregate commodity
% bundles.

% Columns and commodities:
%     01    BREAD     21 POTATOES     41  FURNISH     61 RAILFARE
%     02  CEREALS     22 OTH_VEGS     42 ELEC_APP     62 BUSFARES
%     03 BISCUITS     23    FRUIT     43  OTHHHEQ     63  OTHTRAV
%     04     BEEF     24 OTH_FOOD     44  CONSUMA     64  AUD_VIS
%     05     LAMB     25  CANTEEN     45  PETCARE     65 REC_TOYS
%     06     PORK     26 OTH_SNAC     46  POSTAGE     66 BOOK_NEW
%     07    BACON     27     BEER     47 TELEPHON     67   GARDEN
%     08 POUL_OTH     28 WINESPIR     48 DOMSERVS     68  TVLICEN
%     09     FISH     29     CIGS     49 FEES_SUB     69 ENTERTAI
%     10   BUTTER     30   OTHTOB     50 MENOUTER     
%     11 OIL_FATS     31     RENT     51 WOMOUTER
%     12   CHEESE     32 MORTGAGE     52 KIDOUTER
%     13     EGGS     33 RATESETC     53 OTHCLOTH
%     14 MILKFRES     34   REPAIR     54 FOOTWEAR
%     15 MILKPROD     35      DIY     55 CHEMGOOD
%     16      TEA     36     COAL     56  P_SERVS
%     17   COFFEE     37 ELECTRIC     57 MOTORVEH
%     18 SOFTDRIN     38      GAS     58 MAINTMOT
%     19 SUG_PRES     39  OIL_OTH     59  PET_OIL
%     20 SWE_CHOC     40   FURNIT     60  TAX_INS 

% Assign catagories contingent on the number of classes
if classes == 3
    catagory{1} = 1:26 ;                        %food
    catagory{2} = [27:30 44:45 50:55 64:67] ;   %nondur (including clothing) 
    catagory{3} = [36:39 46:49 56 58:63 68:69]; %services
elseif classes == 4
    catagory{1} = 1:26 ;                        %food
    catagory{2} = [27:30 44:45 55 64:67] ;      %nondur
    catagory{3} = 50:54;                        %clothing
    catagory{4} = [36:39 46:49 56 58:63 68:69]; %services
elseif classes == 5
    catagory{1} = 1:26 ;                        %food
    catagory{2} = [29:30 44:45 55 64:67] ;      %nondur
    catagory{3} = 50:54;                        %clothing
    catagory{4} = [36:39 46:49 56 58:63 68:69]; %services
    catagory{5} = 27:28;                        %alcohol
end
% Columns of non-durables that appear above:
all_nondurables = [1:30 36:39 44:56 58:69]; 
% We exclude columns 31-35, 40-43, 57 as these are durables and not 
% included in any of the catagories.

%% Pre-define variables
p = zeros(periods,classes); % Prices - normalized by median expenditure
wg{periods} = [];           % Shares in period t by household  
totexp{periods} = [];       % Total expenditure in period t by household, 
                            % normalized by median expenditure
Z{periods} = [];

%% Load and clean year-by-year
for t = start_year:end_year
    ss = t-start_year+1;    %Year counter (starts at 1)
  
    % Load data
    data_imp = importdata(['Input/data',  int2str(t), '.csv']);
    data = data_imp.data;
    colheaders = data_imp.colheaders;
    
    % Find columns for cleaning (number cars, kids, ect)
    colcars = find(ismember(colheaders, '"NCARS"')==1);
    colkids = find(ismember(colheaders, '"NUMHHKID"')==1);
    colcouple = find(ismember(colheaders, '"COUPLE"')==1);
    colexp = find(ismember(colheaders, '"TOTEXP"')==1);
    colinc = find(ismember(colheaders, '"HHINCOME"')==1);
 
    % Select couples with one car, more than 0 kids, more than zero income
    % and expenditure
    sel = data(:,colcouple)==1 & data(:,colcars)==1 & data(:,colkids)>0 ...
        & data(:,colexp)>0 & data(:,colinc)>0;
    data = data(sel,:);
    obs = size(data,1);
    
    % Drop families with expenditure or income > 99 percentile
    data = sortrows(data,colinc);
    incomeceiling = data(ceil(0.99.*obs),colinc);
    data = sortrows(data,colexp);
    expceiling = data(ceil(0.99.*obs),colexp); 
    sel = find(data(:,colexp)<=expceiling & data(:,colinc)<=incomeceiling);
    data=data(sel,:);
    
    % Column header for instrument
    colinstrument = find(ismember(colheaders, '"HHINCOME"')==1);
    Z{ss} = data(:,colinstrument);
    
    % Column header numbers for commodities
    colbread = find(ismember(colheaders, '"BREAD"')==1);        %col = 65
    coltob = find(ismember(colheaders, '"OTHTOB"')==1);         %col = 94
    colrent = find(ismember(colheaders, '"RENT"')==1);          %col = 95  
    coldiy = find(ismember(colheaders, '"DIY"')==1);            %col = 99
    colcoal = find(ismember(colheaders, '"COAL"')==1);          %col = 101  
    colentertain = find(ismember(colheaders, '"ENTERTAI"')==1); %col = 134
    
    % Shares: percent of spending on commodity by observation (household)
    % Columns are commodities, rows are observations
    shares = [data(:,colbread:coltob) data(:,colrent:coldiy) data(:,colcoal:colentertain)];
    
    % Define expenditure to be nondurables expenditure
    % w = percent spent on non-durables by household i
    w = sum(shares(:,all_nondurables),2); 
    colexp = find(ismember(colheaders, '"TOTEXP"')==1);
    totexp{ss} = w.*data(:,colexp);
    % Further normalize total expenditure to account for #kids in the
    % household
    kids = data(:,colkids);
    equivscale = 1+(0.29*kids);
    totexp{ss} = totexp{ss}./equivscale;
     
    % Normalize expenditure shares relative to non-durable consumption
    shares(:,all_nondurables) = shares(:,all_nondurables)./kron(ones(1,size(shares(:,all_nondurables),2)),w);
    
    % Aggregate shares into #classes.
    wg{ss} = zeros(size(shares,1),classes);
    for ii = 1:classes
        wg{ss}(:,ii) = sum(shares(:,catagory{ii}),2);
    end
    % Now aggregate shares (wg) sums to one.
 
    % Aggregate Prices (row is commodity and column is year 1975...1999)
    pricefile = 'Input/pall.csv';
    prices_imp = importdata(pricefile);
    prices = prices_imp.data;   

    % One problem is that we do not have good documentation for the price 
    % data, but they state that the ordering is the same as for the 69
    % goods listed above for consumption data

    % Weighting factor: ratio of the average budget share of each category 
    % over the average total budget share for each class of products 
    for ii = 1:classes
        weight = mean(shares(:,catagory{ii}))/sum(mean(shares(:,catagory{ii})));
        % Weighted prices:
        p(ss,ii) = sum(prices(catagory{ii},ss).*(weight'));
    end
    
    % Divide p by 100 (CPI)
    p(ss,:) = p(ss,:)/100;
    
    % Normalize price by median expenditure
    p(ss,:) = p(ss,:)/median(totexp{ss}); 
end

%% Output
clear shares
budgets = p;
shares = wg;
income = totexp;


end











