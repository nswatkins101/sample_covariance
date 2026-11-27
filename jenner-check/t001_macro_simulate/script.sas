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

%simulate(10, 20, 2, 'lognormal');
