/*****************************************************************************
*    Defining and Estimating Effects of Dynamic Treatment Strategies         *
*    Exercise 2 - Estimating Effects                                         *
*    ISPE Course - August 2026                                               *
*    Programmers: Emma McGee & Rienna Russo                                  *
*****************************************************************************/

/* Part 0: Data Setup ------------------------------------------------------ */

/* Load the data */
proc import datafile="pathway_to_your_file\bestmed.csv"
    out=bestmed
    dbms=csv
    replace;
run;


/* Part 1: Inverse Probability Weighting ----------------------------------  */
/* a) ---------------------------------------------------------------------- */
/*   Fit logistic regression models to estimate the denominator              */
/*   probabilities for the non-stabilized IP weights.                        */

/* Weight model for probability of initiating medication at baseline */
proc logistic data=bestmed;
    where Time = 0;
    model strategy(event='1') = 
          /* Baseline time-fixed covariates */
          Female SDI Age Age*Age Insurance
          /* Time-varying covariates */
          eGFR eGFR*eGFR Obesity hba1c
    / link=logit;
    /* Obtain predicted values */
    ods output ParameterEstimates=_parm_tx_b;
    output out=bestmed_base_pred p=p_tx_b;
run;


/* Weight model for probability of continuing medication for SGLT2i strategy */
proc logistic data=bestmed;
    where excused = 0 and strategy = 0 and Time > 0
          and lag1_sglt2i = 1 and lag1_glp1 = 0;
    model sglt2i(event='1') =
          /* Baseline time-fixed covariates */
          Female SDI Age Age*Age Insurance
          /* Time-varying covariates */
          eGFR eGFR*eGFR Obesity hba1c
          /* Time */
          Time Time*Time
    / link=logit;
    /* Obtain predicted values */
    ods output ParameterEstimates=_parm_tx_0;
    output out=bestmed_tx0_pred p=p_tx_0;
run;


/* Weight model for probability of continuing medication for GLP1 strategy */
proc logistic data=bestmed;
    where excused = 0 and strategy = 1 and Time > 0
          and lag1_glp1 = 1 and lag1_sglt2i = 0;
    model glp1(event='1') =
          /* Baseline time-fixed covariates */
          Female SDI Age Age*Age Insurance
          /* Time-varying covariates */
          eGFR eGFR*eGFR Obesity hba1c
          /* Time */
          Time Time*Time
    / link=logit;
    /* Obtain predicted values */
    ods output ParameterEstimates=_parm_tx_1;
    output out=bestmed_tx1_pred p=p_tx_1;
run;

/* Add all predicted values to dataset */
proc sort data=bestmed;          by id Time; run;
proc sort data=bestmed_base_pred;by id Time; run;
proc sort data=bestmed_tx0_pred; by id Time; run;
proc sort data=bestmed_tx1_pred; by id Time; run;

data bestmed_preds;
    merge bestmed          (in=in_base)
          bestmed_base_pred(keep=id Time p_tx_b)
          bestmed_tx0_pred (keep=id Time p_tx_0)
          bestmed_tx1_pred (keep=id Time p_tx_1);
    by id Time;

    if in_base;   /* retain all observations that were in bestmed */

    /* At baseline: p_tx_b only defined at baseline; set missing otherwise */
    if Time ne 0 then p_tx_b = .;

    /* Over follow-up: Set p_tx_0 and p_tx_1 to missing at baseline (Time == 0) */
    if Time = 0 then do;
        p_tx_0 = .;
        p_tx_1 = .;
    end;
run;


/* Calculating time-specific weight contributions in each time period ----- */
data bestmed_wt;
    set bestmed_preds;
    length weight_t 8.;
    if excused = 1 then weight_t = 1;               /* excused */
    else if strategy = 1 and Time = 0 then
        weight_t = 1 / p_tx_b;                      /* initiated GLP1 strategy at baseline */
    else if strategy = 0 and Time = 0 then
        weight_t = 1 / (1 - p_tx_b);                /* initiated SGLT2 strategy at baseline */
    else if strategy = 1 and Time > 0 and glp1 = 1 and sglt2i = 0 and excused = 0 then
        weight_t = 1 / p_tx_1;                      /* following GLP1 strategy, not excused */
    else if strategy = 0 and Time > 0 and sglt2i = 1 and glp1 = 0 and excused = 0 then
        weight_t = 1 / p_tx_0;                      /* following SGLT2 strategy, not excused */
    else weight_t = .; /* weight contributions set to missing at other non-adherent times */
run;


/* c) ---------------------------------------------------------------------- */
/*  Take the cumulative product of the time-specific probabilities estimated */
/*  in (a) and construct time-varying IP weights for each eligible individual*/

/* Make sure rows are ordered from first to last time within ID */
proc sort data=bestmed_wt;
    by id Time;
run;

/* Take the cumulative product of the time-specific weight contributions */
/* Beginning from the first time (start of follow-up), separately for each ID */
data bestmed_wt;
    set bestmed_wt;
    by id;
    retain w_a;
    if first.id then w_a = .;

    if weight_t ne . then do;
        if first.id then w_a = weight_t;
        else w_a = w_a * weight_t;
    end;
    else w_a = .;
run;


/* Part 2: Censoring ------------------------------------------------------- */
/* b) ---------------------------------------------------------------------- */
/*   Implement censoring rules for each treatment strategy to create dataset */
/*   in which all eligible individuals follow 'assigned' strategy.           */

/* Create variable to censor people */
data bestmed_dev;
    set bestmed_wt;
    if (   (strategy = 1 and glp1 ne 1 and excused = 0)
        or (strategy = 1 and glp1 = 1 and sglt2i = 1 and excused = 0)
        or (strategy = 0 and sglt2i ne 1 and excused = 0)
        or (strategy = 0 and sglt2i = 1 and glp1 = 1 and excused = 0) ) then
        dev_time2 = Time;
    else dev_time2 = 999; /* dev_time2 is in intermediary variable */
run;
proc summary data=bestmed_dev nway;
    class id;
    var dev_time2;
    output out=cens_time_dynamic(drop=_type_ _freq_)
           min=cens_time_dynamic; /* take the minimum of dev_time2 within id, i.e., the first time they deviated */
run;

/* Truncate follow-up once they stopped following their 'assigned' strategy */
proc sort data=bestmed_wt;        by id Time; run;
proc sort data=cens_time_dynamic; by id;      run;

data bestmed_cens;
    merge bestmed_wt(in=a) cens_time_dynamic(in=b);
    by id;
    if a;
    if Time < cens_time_dynamic;  /* truncate follow-up */
run;

/* Check censoring times */
proc freq data=bestmed_cens;
    where Time = 0;
    tables cens_time_dynamic / nocum nopercent;
run;

/* Check weights and censored dataset construction */
proc print data=bestmed_cens (obs=100);
    var id Time cens_time_dynamic strategy excused
        glp1 sglt2i p_tx_b p_tx_1 p_tx_0 weight_t w_a;
run;


/* c) ---------------------------------------------------------------------- */
/*   Truncate IP weights at the 99th percentile.                             */

/* Summarize weight distribution */
proc univariate data=bestmed_cens;
    var w_a;
    output out=wquant pctlpts=99 pctlpre=p_;
run;

/* Truncate weights */
data _null_;
    set wquant;
    call symputx('p99_w', p_99);
run;

data bestmed_cens;
    set bestmed_cens;
    if w_a >= &p99_w then w_a = &p99_w;
run;

/* Summarize weight distribution after truncating */
proc univariate data=bestmed_cens;
    var w_a;
run;


/* Part 3: Dynamic Marginal Structural Model ------------------------------- */
/* a) ---------------------------------------------------------------------- */
/*   Fit an IP weighted pooled logistic regression model for the outcome of  */
/*   MACE                                                                    */

data bestmed_cens;
    set bestmed_cens;
    Time2 = Time*Time;
run;

proc logistic data=bestmed_cens outmodel=msm_model;
    model mace(event='1') = Time strategy Time2
                            Time*strategy Time2*strategy
    / link=logit;
    weight w_a;
run;


/* b) ---------------------------------------------------------------------- */
/*   Using the parameter estimates from the model in (a), estimate the       */
/*   marginal risk of MACE at each month of follow-up under each strategy.   */

/* Create dataset of all timepoints (0-47) under each strategy (0,1) */
data risk_results;
    do strategy = 0, 1;
        do Time = 0 to 47;
            Time2 = Time*Time;
            output;
        end;
    end;
run;

/* # Estimate discrete-time hazards under each strategy */
proc logistic inmodel=msm_model;
    score data=risk_results out=risk_results_scored(rename=(P_1=hazard));
run;

/* Estimate survival from cumulative product of (1-hazard) for each strategy */
/* Estimate risks from survival probabilities */
proc sort data=risk_results_scored;
    by strategy Time;
run;

data risk_results_scored;
    set risk_results_scored;
    by strategy;
    retain survival;
    if first.strategy then survival = 1;
    survival = survival * (1 - hazard);
    risk = 1 - survival;
run;


/* c) ---------------------------------------------------------------------- */
/*   Using risks from (b), construct adjusted risk curves for the outcome of */
/*   MACE under each strategy.                                               */

/* Prepare data */
/* Shift Time by +1 (outcomes appearing in interval k represent those that happened in interval k+1) */
data risk_plot;
    set risk_results_scored;
    Time = Time + 1;
run;

/* Add Time=0, risk=0 for each strategy */
data risk_plot;
    set risk_plot
        end=eof;
    output;
run;

data risk_plot;
    set risk_plot
        (where=(1=1))
        end=eof;
    output;
    if eof then do;
        Time = 0; risk = 0;
        strategy = 0; output;
        strategy = 1; output;
    end;
run;

proc sort data=risk_plot;
    by strategy Time;
run;

/* Plot risk curves */
proc format;
    value stratf
        0 = 'SGLT2i' 
        1 = 'GLP-1RA'; /* Labels for legend */
run;

proc sgplot data=risk_plot;
    format strategy stratf.;             /* Use labels from above */
    styleattrs datacontrastcolors=(cx000080 cx56B4E9) ;

    series x=Time y=risk /
        group=strategy
        grouporder=ascending
        name="risk";

    xaxis label="Months" values=(0 to 48 by 6);
    yaxis label="Risk (%)"
          min=0 max=0.075
          values=(0 0.025 0.05 0.075)
          valuesdisplay=("0.0%" "2.5%" "5.0%" "7.5%");

    keylegend "risk" /
        title="Strategy"
        location=outside
        position=bottom
        valueattrs=(family=Arial);
run;


/* d) ---------------------------------------------------------------------- */
/*   Using the risks from (b), estimate the 4-year risk difference           */
/*   and ratio for MACE                                                      */

/* # 4-year risks, risk difference, and risk ratio  */
proc sql noprint;
    select risk into :risk1
    from risk_results_scored
    where Time = 47 and strategy = 1;

    select risk into :risk0
    from risk_results_scored
    where Time = 47 and strategy = 0;
quit;

data results;
    risk1 = &risk1;
    risk0 = &risk0;
    rd    = risk1 - risk0;        /* risk difference */
    rr    = risk1 / risk0;        /* risk ratio      */
run;

proc print data=results;
    format risk1 risk0 rd rr 8.4;
run;
