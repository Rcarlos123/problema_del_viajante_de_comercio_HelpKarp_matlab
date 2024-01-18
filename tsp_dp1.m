
    function [OptimalTour,mincost]=tsp_dp1(cities, Dmatrix)

    % Initialization Process
    if nargin==0
        cities=[0,2,9,10;1,0,6,4;15,7,0,8;6,3,12,0] 

       % cities=random('unif',0,100,10,2);
    end
    [NumOfCities,dummy]=size(cities);
    %generador de numeros primos
    Primes=primes(NumOfCities*10);
    
    %Primes=[2 3 5 7 11 13 17 19 23 29 31 47]

    if nargin<2 %No Dmatrix used 
        %se crea la matriz 
        D=diag(inf*ones(1,NumOfCities)); %Asigne un costo infinito a viajar de una ciudad a sí misma
        for i=1:NumOfCities %Asignar las distancias entre pares de ciudades
            for j=i+1:NumOfCities
                D(i,j)=norm(cities(i,:)-cities(j,:));
                D(j,i)=D(i,j);
            end
        end
    else
        D=Dmatrix;
    end
    NumOfDataSets=1;
   
    for i=2:NumOfCities
        NumOfDataSets=NumOfDataSets+nchoosek(NumOfCities,i);
    end
    
    
    Data(NumOfDataSets).S=[];
    Data(NumOfDataSets).l=0;
    Data(NumOfDataSets).cost=inf;
    Data(NumOfDataSets).pre=[];
    Data(NumOfDataSets).m=[];
    LookUpTable(NumOfDataSets)=0;
    %Definir una estructura de datos que contenga los siguientes datos que necesitamos
    %Para luego. Esta estructura de datos utiliza la misma notación utilizada en el documento.
    % por Held y Karp (1962):
    % S - el conjunto de ciudades en el recorrido.
    % l - la última ciudad visitada en el conjunto S.
    % cost - el costo de un recorrido, que incluye toda la ciudad en S y termina en l.
    %Además, los siguientes elementos de datos se utilizan en el conjunto de datos para reducir
    % tiempo de ejecución:
    % Pre - el índice del conjunto de datos predecesor, es decir, el que tiene Spre=S-{l}
    % m - la ciudad en S-{l} que produjo el costo más bajo C(Spre,m)+D(m,l).
    % Este índice facilitará la generación del tour óptimo sin
    % cálculos adicionales.
    Data(1).S=[1];
    Data(1).l=1;
    Data(1).cost=0;
    Data(1).Pre=[];
    Data(1).m=[];
 
    for s=2:NumOfCities
        Data(s).S=[Data(1).S,s];
        Data(s).l=s;
        Data(s).cost=D(s,1);
        Data(s).Pre=1;
        Data(s).m=1;
        LUT=calcLUT(Data(s).S,s,Primes);
        LookUpTable(s)=LUT;
    end
    IndexStartPrevStep=2;%indice en Datos que marcan el comienzo del paso anterior
    IndexLastStep=NumOfCities; %indice en Datos que marcan el final del paso anterior
    CurrentData=IndexLastStep; %indice en Datos que marcan el conjunto de datos actual

    %Este es el bucle principal de programación dinámica
    for s=3:NumOfCities
       
        TempSets=nchoosek(2:NumOfCities,s-1);

       
        NumOfSets=size(TempSets);

        for j=1:NumOfSets(1)
           
            for k=1:NumOfSets(2)

                SminuskSet=[1,TempSets(j,1:k-1),TempSets(j,k+1:NumOfSets(2))];%este es el conjunto S-{k} 
                
                
                candidatecost(2:length(SminuskSet))=inf;
                
               
                indices=[];
                for mm=2:length(SminuskSet) % elige una ciudad en S-{k} que será la última
                   
                    LUV=calcLUT(SminuskSet,SminuskSet(mm),Primes);
                  
                    index=find(LUV==LookUpTable(IndexStartPrevStep:IndexLastStep));
                    
                    index=index+IndexStartPrevStep-1;
                   
                    if index==0
                        candidatecost(mm)=inf
                        
                    else
                        candidatecost(mm)=Data(index).cost+D(SminuskSet(mm),TempSets(j,k));
                      
                        indices(mm)=index;
                    end
                  
                end
                [mincost,indexcost]=min(candidatecost(2:end));
                CurrentData=CurrentData+1;
                Data(CurrentData).S=[1,TempSets(j,:)];
                Data(CurrentData).l=TempSets(j,k);
                Data(CurrentData).cost=mincost;
                Data(CurrentData).Pre=indices(indexcost+1);
                Data(CurrentData).m=SminuskSet(indexcost+1);
                LookUpTable(CurrentData)=calcLUT(Data(CurrentData).S,TempSets(j,k),Primes);

            end
       
        end
        IndexStartPrevStep=IndexLastStep+1;
        IndexLastStep=CurrentData;    

    end
    mm=0;
    % Ahora agregue la distancia desde la última ciudad hasta la ciudad 1
    for i=IndexStartPrevStep:IndexLastStep
        mm=mm+1;
        candidatecost(mm)=Data(i).cost+D(Data(i).l,1);
    end
    %Encuentre el que minimice la distancia total
    [mincost,indexcost]=min(candidatecost);
    Temp=Data(IndexStartPrevStep+indexcost-1);
    % Genere el recorrido óptimo al regresar de la última ciudad a su
    %antecesores
    OptimalTour=1;
    while ~isempty(Temp.Pre)
        OptimalTour=[OptimalTour,Temp.l];
        Temp=Data(Temp.Pre);
    end

    OptimalTour=[OptimalTour,1];


    function LUT=calcLUT(vec,last,Primes)
   
    j=length(vec);
  
    LUT=Primes(last);
    

    for i=2:j
        LUT=LUT*Primes(vec(i));
       

    end
