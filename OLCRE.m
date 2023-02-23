function OLCRE
    clear,clc;
    close all;
   
    %Input network bilateral information
    %Take the karate network for example
    str = './dataset/karate/karate_2Link.dat'; 
       
    global linkData2;
    linkData2 = load(str);    
    
    w=sparse(linkData2(:,1),linkData2(:,2),ones(size(linkData2,1),1));
    G=graph(w);
  
    count = max(max(linkData2(:,1)),max(linkData2(:,2)));
    
    global nodesdegree;   
    nodesdegree = zeros(count,2);
    for i=1:count
        nodesdegree(i,1) = i;
        temp = find(linkData2(:,2)==i);
        nodesdegree(i,2) = length(temp);
    end
    
    nodes = nodesdegree(:,1);
       
    %Select seed communities
    LCR = zeros(count,3); 
    for i =1:length(LCR)
        LCR(i,1)=i;
    
        n = [];
        n=neighbors(G,LCR(i,1));   
        
        P2 = size(n,1);      
        for z=1:length(n)
            P2 = P2 + size(neighbors(G,n(z)),1);        
        end
        
        sum = 0;
        sumI = 0;
        for j=1:length(n)  
            sum = sum + size(intersect(n,neighbors(G,n(j))),1);      
        end
        sumI = sum / 2 + size(n,1);        
        
        sumE = 0;
        n = union(n,LCR(i)); 
        for j=1:length(n)
            sumE = sumE + size(setdiff(neighbors(G,n(j)),n),1);         
        end
        P1 = sumI / sumE;    
        
        LCR(i,2) = P1*P2;      
     
    end 
    ULCR=[]; 
    for i = 1:length(LCR)
        if LCR(i,3)==1      
            continue;
        end 
        LCR(i,3)=1;
        linJu = [];
        linJu=neighbors(G,LCR(i,1));
        flag = 0;   
        value = LCR(i,1);   
        uValue = LCR(i,1);  
         
        while true
                for j =1 :length(linJu)
                    LCR(linJu(j),3)=1;
                    if LCR(linJu(j),2)>LCR(value,2) 
                            value = LCR(linJu(j));  
                            continue;
                    else
                            if j == length(linJu)      

                            else 
                                    continue;

                            end    
                        
                    end    
                    
                    if value == uValue 
                        flag = 1;
                    else    
                        linJu = neighbors(G,LCR(value,1));  
                        uValue = value;
                    end    
                end
                if flag ==1;       
                    ULCR = union(ULCR,LCR(uValue));    
                    break;
                end  
        end  
    end   
    
    ULCR = ULCR';
 
    CN = zeros(length(ULCR),count);
    DN = [];
    MN=[];
    ON=[];
   
    %Seed communities for expansion
    for i=1:length(ULCR)     
        seedN = [];    
        neighN = [];   
        bN = [];       
        tcg = 0;      
        sumW = 0;   
        
        seedN = [neighbors(G,ULCR(i,1));ULCR(i,1)];     
        seedN = sortrows(seedN);  
    
        neighN = [];
        for h=1:size(seedN,1)
            nnN = setdiff(neighbors(G,seedN(h)),seedN);    
            neighN = union(neighN,nnN);    
        end
        neighN;       
        bN = [];
        for j=1:size(seedN,1)
            if isempty(setdiff(neighbors(G,seedN(j)),seedN))==0     
                bN = [bN;seedN(j)];
            end   
        end
        bN;    
        neigh = [];   
        for k=1:length(seedN)   
            neigh = intersect(neighbors(G,seedN(k,1)),seedN);
            sum = 0;
            for q=1:length(neigh)
                sum = sum + size(intersect(neigh,neighbors(G,neigh(q))),1);    
            end
            tcg = tcg + sum / 2 ;       
        end
        tcg = tcg / 3; 
        for l = 1:length(seedN)
            sumW = sumW + length(setdiff(neighbors(G,seedN(l)),seedN)); 
        end
        sumW;      
        x = 1;
        while true
            if isempty(neighN)==1   
                break;    
            end     
            tcg2 = 0;   
            nh = [];   
            nh = intersect(seedN,neighbors(G,neighN(x)));
            sum2 = 0;
            for u=1:length(nh)
                sum2 = sum2 + size(intersect(nh,neighbors(G,nh(u))),1);         
            end
            tcg2 = tcg + sum2 / 2 ;   
            sumW2 = 0;  
            sumW2 = sumW - length(intersect(neighbors(G,neighN(x)),seedN)) + length(setdiff(neighbors(G,neighN(x)),seedN)); 
            I = tcg / length(seedN);
            E = sumW / length(bN);
            I2 = tcg2 / (length(seedN) + 1);
            if length(setdiff(neighbors(G,neighN(x)),seedN))==0 
                E2 = sumW2 / (length(bN) - 1);
            else  
                E2 = sumW2 / length(bN);
            end   
            if I2 - I > 0 && E2 - E < 0
                seedN = union(seedN,neighN(x));   
                nnN = [];
                nnN = setdiff(neighbors(G,neighN(x)),seedN);     
                neighN = union(setdiff(neighN,neighN(x,1)),nnN);  
                tcg = tcg2;
                sumW = sumW2;
                x = 1;
            else
                x = x + 1;
                if x>length(neighN)
                    break;
                end    
            end   
        end  
        for z = 1:length(seedN)
            CN(i,z) = seedN(z);  
        end  
        for v = 1:size(CN,1)
           DN = union(DN,CN(v,:));  
        end
        DN = [DN;zeros(count-length(DN),1)]; 
    end 
    
    

    
    
    %Handle missing nodes
    
    if length(setdiff([nodes;0],DN)) == 0 
        fprintf('All nodes have been traversed\n');
    else
        
        MN = setdiff([nodes;0],DN);     
    end 

    K = zeros(size(CN,1),2);
    for i =1:length(MN)    
        for j =1 :size(CN,1)
            sumNB1 = 0;
            sumWB1 = 0;
            sumNB2 = 0;
            sumWB2 = 0;
            
            U = CN(j,:);
            U(find(U==0)) = [];
            
            for n =1:length(U)
                sumNB1 = sumNB1 + length(intersect(neighbors(G,U(n)),U));
                sumWB1 = sumWB1 + length(setdiff(neighbors(G,U(n)),U));
            end
            sumNB1 = sumNB1 / 2;
            
            U = union(U,MN(i));
            
            for m =1:length(U)
                sumNB2 = sumNB2 + length(intersect(neighbors(G,U(m)),U));
                sumWB2 = sumWB2 + length(setdiff(neighbors(G,U(m)),U));
            end
            sumNB2 = sumNB2 / 2;
            
            K(j,1) = j;
            K(j,2) = (sumNB2 / sumWB2) - (sumNB1 / sumWB1);   
        end 
        
        [row,col] = find(K==max(K(:,2))); 
        [row2,col2 ]= find(CN(row,:)==0);  
        CN(row,col2(1)) = MN(i);  
          
    end  
    


    
    %Simplify and optimize overlapping nodes
    for i=1:(size(CN,1)-1)
        for j = (i+1):size(CN,1)
        ON = union(ON,intersect(CN(i,:),CN(j,:)));
        ON = setdiff(ON,0);
        end
    end     

    for i = 1:size(CN,1)
       for j = 1:length(ON) 
           if intersect(CN(i,:),ON(j)) == ON(j)
               
               sumNB1 = 0;
               sumWB1 = 0;
               sumNB2 = 0;
               sumWB2 = 0;
               U=[]; 
               U = setdiff(CN(i,:),ON);    
               U(find(U==0)) = []; 
               
               for n =1:length(U)
                    sumNB1 = sumNB1 + length(intersect(neighbors(G,U(n)),U));
                    sumWB1 = sumWB1 + length(setdiff(neighbors(G,U(n)),U));
               end
               sumNB1 = sumNB1 / 2;
               
               U = union(U,ON(j));
               for m = 1:length(U)
                   sumNB2 = sumNB2 + length(intersect(neighbors(G,U(m)),U));
                   sumWB2 = sumWB2 + length(setdiff(neighbors(G,U(m)),U));
               end
               sumNB2 = sumNB2 / 2;
               
               if (sumNB2 / sumWB2) - (sumNB1 / sumWB1) < 0 
                   CN(i,find(CN(i,:)==ON(j)))=0;                
               end    
           end    
       end          
    end 

    for i =1:size(CN,1)
        C=[];
        C = union(C,CN(i,:));
        C = setdiff(C,0);
        C = C';
        CN(i,:) = [C,zeros(1,(count - length(C)))];
        
    end    
    
   %When the M values of the overlapping node corresponding to the communities where it is located are all negative, it will be added to the community with the corresponding largest M value. 
    UMN = []; 
    for i=1:size(CN,1)
        UMN = union(UMN,CN(i,:));
        UMN = setdiff(UMN,0);
    end
    UMN = setdiff(nodes,UMN);
    
    if isempty(UMN)==0
        KK = zeros(size(CN,1),2); 

        for i =1:length(UMN)    
            for j =1 :size(CN,1) 
                sumNB1 = 0;
                sumWB1 = 0;
                sumNB2 = 0;
                sumWB2 = 0;

                U = CN(j,:);
                U(find(U==0)) = []; 

                for n =1:length(U)
                    sumNB1 = sumNB1 + length(intersect(neighbors(G,U(n)),U));
                    sumWB1 = sumWB1 + length(setdiff(neighbors(G,U(n)),U));
                end
                sumNB1 = sumNB1 / 2;

                U = union(U,UMN(i));

                for m =1:length(U)
                    sumNB2 = sumNB2 + length(intersect(neighbors(G,U(m)),U));
                    sumWB2 = sumWB2 + length(setdiff(neighbors(G,U(m)),U));
                end
                sumNB2 = sumNB2 / 2;

                KK(j,1) = j;
                KK(j,2) = (sumNB2 / sumWB2) - (sumNB1 / sumWB1)
            end 
            [row,col] = find(KK==max(KK(:,2)));  
            [row2,col2 ]= find(CN(row,:)==0);  
            CN(row,col2(1)) = UMN(i);   
        end  
        
    end    
    
    fprintf('Save the result to a file\n');
    doc='./dataset/karate/karate_OLCRE';
    saveTxt(CN,doc);
    

    
end 





%A function that saves community results to a file
function saveTxt(labels,name) 
    str = [name '.dat'];
	fid = fopen(str,'w');
	for i = 1 : size(labels,1)
        for j = 1 : size(labels,2)
            fprintf(fid,'%d ',labels(i,j));
        end
        fprintf(fid,'\n');
    end
	fclose(fid); 
	fprintf('Create file %s is OK!\n',str);
end