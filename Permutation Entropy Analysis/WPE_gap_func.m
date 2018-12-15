function [ H, n ] = WPE_gap_func( x, t, wl, tau )
%WPE: Garland weighted permutation entropy 
%   Estimates the entropy of any real-valued time series
%   Excludes permutations covering time gap.
%   Excludes permutations covering a string of zeros.

% %To troubleshoot code:
%     x=(CHA+CHJ);
%     t=mtime;
%     tau=1;
%     wl=3;

xl=x;
    for m=1:wl-1
        xl=[xl vertcat(x(tau*m+1:end),zeros(tau*m,1))]; %table of lags
    end
    xl(end-((wl-1)*tau):end,:)=[];% remove rows with zeros at the end of xl 
    
    %Remove permutations covering gaps in the time series
        tgap=vertcat(0,diff(t)>=2);%find data gaps
        startremove=find(tgap==1)-wl;
        endremove=find(tgap==1);
        removeindex=zeros(length(startremove),wl+1);
        for r=1:length(startremove)
            removeindex(r,:)=startremove(r):1:endremove(r);
        end
        ri=unique(sort(reshape(removeindex,(wl+1)*length(removeindex),1)));
        ri(ri<=0)=[]; %remove zeros and negative values
        ri(ri>length(xl))=[];
        xl(ri,:)=[]; 
         n=length(xl);      
        
            xltemp=xl(:,1:wl); %just the columns of xl for this word length
            maxH=sum(log2(2:wl)); %maximum entropy
            [~,ind]=sort(xltemp,2);
            word=ind*(10.^(0:wl-1)'); %list of words in the time series
            vword=var(xltemp,[],2); %calculate the weight of each word as the variance of the word elements
            u=unique(word);
            for k=1:length(u) %for each unique word
                inc=(word==u(k)); %index the matching words in the time series word list
                vw(k)=sum(inc.*vword); %sum the variances of that word
            end
          
            W=vw/sum(vw); %the weighted probability of a permutation
                       
             %if any of the variances are zero then this returns W=NaN when it should be 1
            if vw==0
                W=1;
            end
            
            Hobs=sum(W.*log2(W)); %the weighted permutation entropy
            
            H=-Hobs/maxH; %normalized WPE relative to the maximum entropy
            clear vw k inc
end



