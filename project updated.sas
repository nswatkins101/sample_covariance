/*
Final Project
Wishart Distribution, Marchenko-Pastur Law
Nathaniel Watkins
Jing Chen
*/


options mprint mlogic;

*** DEFINE ALL MACROS;

* a macro that generates samples of correlated random variables, with some fixed covariance structure
* which is hard-coded. This code is just used to generate an example of correlations.;
%macro simulate_correlated(num_samples, sample_size, dist)
	* generate a sample of size n in 2 variables;
	data samples (drop=d);
		do sampleid=1 to &num_samples;
			do n=1 to &sample_size;
				do d=1 to 2;
					* generate a random variable from the specified distribution;
					x=rand(&dist);
					output;
				end;
			end;
		end;
	run;

	proc transpose data=samples out=samples_trans (drop=_NAME_) prefix=x;
		by sampleid n;
	run;

	* plot sample using 2d scatterplot;

	proc sgplot data=samples_trans (where=(sampleid=1));
		scatter y=x2 x=x1;
		title "Uncorrelated Sample with n=&sample_size from &dist distribution";
	run;

	* apply a linear transform to get correlated random variables;

	data corr_data (drop=x1 x2);
		set samples_trans;
		y1=(3*x1+x2)/sqrt(2);
		y2=(3*x1-x2)/sqrt(2);
	run;

	proc sgplot data=corr_data (where=(sampleid=1));
		scatter y=y2 x=y1;
		title "Correlated Sample with n=&sample_size. from &dist. distribution";
	run;
%mend;

* a macro to generate the simulation data such that the covariance structure is as given;
%macro simulate(num_samples, sample_size, dim, dist);
	* generate a sample of size sample_size in d variables;
	* calculate the standard deviation of the distribution;
	data test_sample;
		do n=1 to 10000;
			* generate a random variable from the distribution;
			x=rand(&dist.);
			output;
		end;
	run;
	proc sql;
		select std(x)
			into :std_dev
		from test_sample;
	quit;
	
	
	data samples (drop=d);
		do sampleid=1 to &num_samples;
			do n=1 to &sample_size;
				do d=1 to &dim;
					* generate a random variable from the distribution;
					x=rand(&dist.)/&std_dev.;
					output;
				end;
			end;
		end;
	run;

	proc transpose data=samples out=samples_trans (drop=_NAME_) prefix=x;
		by sampleid n;
	run;

	* plot sample using 2d scatterplot;

	proc sgplot data=samples_trans (where=(sampleid=1));
		scatter y=x2 x=x1;
		title "Sample with n=&sample_size from &dist distribution";
	run;

	* rename dataset;
	data corr_data (drop=x:);
		set samples_trans;
		%do i=1 %to &dim;
			y&i. = x&i.;
		%end;
	run;

%mend simulate;

* a macro to calculate and plot the empirical spectral measure;
%macro esm(num_samples, sample_size, dim, dist);
	* calculate covariance matrix of sample;

	proc corr data=corr_data outp=OutCorr noprint nomiss cov;
		var y:;
		by sampleid;
	run;

	* use SAS/IML;

	proc iml;
		use OutCorr where(_TYPE_="COV");
		do i=1 to &num_samples;
			read all var("y1":"y&dim.") into cov[colname=varNames] where(sampleid=i);
			* calculate eigenvalues;
			sample_val=eigval(cov);
			val=val // sample_val;
		end;
		* write total_val to a SAS dataset;
		create eigenvalues var {val};
		append;
		close eigenvalues;
	quit;
	* create histogram from eigenvalues;
	title "Histogram of eigenvalues n=&sample_size. from &dist. distribution";
	ods graphics on;

	proc univariate data=eigenvalues (rename=(val=eigenvalue)) noprint;
		where eigenvalue between 0 and 5;
		histogram eigenvalue / nmidpoints=100

		/*rtinclude*/
		outhistogram=hist odstitle=title;
	run;

	* create spectral density from histogram data, sgplot has more flexibility than proc univariate histogram;

	/*data spectral_density (drop = _: binsize);
	set hist;
	binsize = 0.1;
	eigenvalue = _MIDPT_;
	density = _count_/binsize;
	run;
	
	proc sgplot data=spectral_density;
	scatter y=density x=eigenvalue;
	run;*/

%mend esm;

* now create macro that creates report, will loop over sample size and several distributions
  to generate an interesting report.;
  
%macro sample_size_exp(dist);

	%do samplesize= 10 %to 510 %by 100;
		%simulate(10000, &samplesize, 2, &dist);
		%esm(10000, &samplesize, 2, &dist);
	%end;

%mend;

* a different type of experiment would increment the dimensionality of the problem as well
  so that the ratio of dimensionality to sample_size is fixed, <1. This results in the
  Marchenko-Pastur Law.;
  
%macro mp_exp(ratio, dist);
	%do sample_size= 20 %to 500 %by 20;
		* ceiling so that dimensionality is integer and non-zero;
		%let dim=%sysfunc(ceil(&ratio.*&sample_size.));
		%simulate(1000, &sample_size., &dim., &dist.);
		%esm(1000, &sample_size., &dim., &dist.);
	%end;

%mend;

* the above just shows the Marchenko-Pastur Law is the asymptotic spectral density, we
  now want to show that it is universal. We plot the spectral densities together in one plot,
  they should approach each other.;
  
  
  
* Jing's code;

%macro universality(dist1=normal, dist2=lognormal, num_samples=1000, sample_size=100, ratio=0.5);
   %let dim = %sysfunc(ceil(&ratio.*&sample_size.));  /* Calculate the dimensionality for the ratio=d/n */
   %put &dim.;
   /* Generate samples and calculate the eigenvalues for the first distribution */
   %simulate(&num_samples., &sample_size., &dim., &dist1.);
   %esm(&num_samples., &sample_size., &dim., &dist1.);
   data eigenvalues1(rename=(val=eigenvalue1));
       set eigenvalues;
   run;
   /* Generate samples and calculate the eigenvalues for the second distribution */
   %simulate(&num_samples., &sample_size., &dim., &dist2.);
   %esm(&num_samples., &sample_size., &dim., &dist2.);
   data eigenvalues2(rename=(val=eigenvalue2));
       set eigenvalues;
   run;
   /* Merge the two sets of eigenvalues */
   data combined;
       merge eigenvalues1 eigenvalues2;
   run;
   /* Plot the spectral densities together */
   proc sgplot data=combined;
       where eigenvalue1 between 0 and 5;
       where eigenvalue2 between 0 and 5;
       histogram eigenvalue1 / binwidth=0.1 transparency=0.5;
       histogram eigenvalue2 / binwidth=0.1 transparency=0.5;
       density eigenvalue1 / type=kernel;
       density eigenvalue2 / type=kernel;
   run;
%mend universality;

*** TEST MACROS;
%simulate(10, 20, 2, 'lognormal');

%sample_size_exp('normal');
%sample_size_exp('lognormal');

%mp_exp(0.1, 'normal');
%mp_exp(0.1, 'lognormal');
/* Call the macro for the universality */
%universality(dist1='normal', dist2='lognormal', num_samples=1000, sample_size=20, ratio=0.1);
%universality(dist1='normal', dist2='lognormal', num_samples=100, sample_size=2000, ratio=0.1);












