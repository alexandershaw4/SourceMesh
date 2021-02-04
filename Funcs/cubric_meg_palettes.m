function [map2,map3,map4] = cubric_meg_palettes
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:255,
	  map2(i,1)=0;
	  map2(i,2)=0;
	  map2(i,3)=0;
	  map4(i,1)=0;
	  if i>=128,
		map4(i,1)=0.5+0.5*((i-128)/127.0)^0.6;
	  	map2(i,1)=3*(i-128)/128.0;
	  	map2(i,2)=3*0.7*(i-128)/128.0;
	  end
	  if i<128,
		map4(i,1)=0.5-0.5*((127-i)/126.0)^0.6;
	  	map2(i,2)=3*0.7*(128-i)/128.0;
	  	map2(i,3)=3*(128-i)/127.0;
	  end
	  if map4(i,1)>1.0, map4(i,1)=1.0; end
	  if map4(i,1)<0.0, map4(i,1)=0.0; end

	  if map2(i,1)>1.0, map2(i,1)=1.0; end
	  if map2(i,2)>0.7, map2(i,2)=0.7; end
	  if map2(i,3)>1.0, map2(i,3)=1.0; end
	  map4(i,2)=map4(i,1);
	  map4(i,3)=map4(i,1); 	
	end
	M=1;
	for i=128:255,
			if i>=128&&i<=149,
				r=M;g=M-0.0476*M*(i-128);b=M;
			end			
			if(i>=150&&i<=170),
				r=M-0.05*M*(i-150);g=0;b=M;
			end
			if(i>=171&&i<=191),
				r=0;g=0;b=M-0.05*M*(i-171);
			end
			if(i>=192&&i<=212),
				r=0.05*M*(i-192);g=0;b=0;
			end
			if(i>=213&&i<=232),
				r=M;g=0.029*M*(i-213);b=0.03*M;
			end
			if(i>=233&&i<=255),
				r=M;g=0.59*M+0.0195*M*(i-233);b=0.03*M;
			end
		if r<0,
		  r=0;
		end
		if g<0,
		  g=0;
		end
		if b<0,
		  b=0;
		end
		if r>1,
		  r=1;
		end
		if g>1,
		  g=1;
		end
		if b>1,
		  b=1;
		end
		map3(i-127,1)=r;
		map3(i-127,2)=g;
		map3(i-127,3)=b;
    end
end

