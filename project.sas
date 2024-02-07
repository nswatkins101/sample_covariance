/*
Final Project
Wishart Distribution, Marchenko-Pastur Law
*/

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

* a macro to generate the simulation data such that the covariance structure identity;
%macro simulate(num_samples, sample_size, dim, dist);
	* generate a sample of size n in d variables;

	data samples (drop=d);
		do sampleid=1 to &num_samples;

			do n=1 to &sample_size;

				do d=1 to &dim;
					* generate a random variable from the normal distribution;
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
		title "Sample with n=&sample_size from &dist distribution";
	run;
	
	* rename the dataset;

%mend simulate;

* a macro to calculate and plot the empirical spectral measure;
%macro esm(num_samples, sample_size, dim, dist);
	* calculate covariance matrix of sample;

	proc corr data=corr_data outp=OutCorr noprint nomiss cov;
		var y1 y2;
		by sampleid;
	run;

	* use SAS/IML;

	proc iml;
		do i=1 to &num_samples;
			use OutCorr where(sampleid=i & _TYPE_="COV");
			read all var {y1 y2} into cov[colname=varNames];
			print cov;
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
		title "Histogram of Eigenvalues n=&sample_size. from &dist. distribution";
		ods graphics on;

	proc univariate data=eigenvalues (rename=(val=eigenvalue)) noprint;
		where eigenvalue between 0 and 25;
		histogram eigenvalue / nmidpoints=100

		/*rtinclude*/
		outhistogram=hist odstitle=title;
	run;

	* create spectral density from histogram data, code not used;

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

*%esm(1000, 25, 2, 'normal');

* now use non-trivial mv dist, there is still correlation because of the linear transform;
*%esm(1000, 25, 2, 'lognormal');

* I will use a different univariate distribution and induce correlation by the same linear transform;
* calculate the sample covariance of this sample;
* generate the empirical spectral distribution of the sample covariance;
* calculate the eigenvalues of the sample cov and write these eigenvalues to a data set;
* plot a histogram of the eigenvalues, the shape should be like Wishart;

* now create macro that creates report, will loop over sample size and several distributions
  to generate an interesting report.;
  
%macro sample_size_exp(dist);

	%do samplesize= 80 %to 90 %by 10;
		%simulate(1000, &samplesize, 2, &dist);
		%esm(1000, &samplesize, 2, &dist);
	%end;

%mend;

%sample_size_exp('normal');
%sample_size_exp('lognormal');

* a different type of experiment would increment the dimensionality of the problem as well
  so that the ratio of dimensionality to sample_size is fixed, <1. This results in the
  Marchenko-Pastur Law.;
  
%macro mp_exp(ratio, dist);
	
	%do sample_size= 100 %to 500 %by 100;
		* ceiling so that dimensionality is integer and non-zero;
		%simulate(1000, sample_size, ceiling(ratio*sample_size), &dist);
		%esm(1000);
	%end;

%mend;

%mp_exp(0.1, 'normal');

* the above just shows the Marchenko-Pastur Law is the asymptotic spectral density, we
  now want to show that it is universal. We plot the spectral densities together in one plot,
  they should approach each other.;
  
  
  


