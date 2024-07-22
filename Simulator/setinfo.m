function setinfo(~,~)

global data 

persistent ct
if isempty(ct)
    for i = 1:data.numUEs
        ct(i)  = 1;
    end
end

persistent rxvalSave
if isempty(rxvalSave)
    for i = 1:data.numUEs
        rxvalSave(i)  = 0;
    end
end

persistent txvalSave
if isempty(txvalSave)
    txvalSave  = 0;
end


rx = zeros(1,data.numUEs);
for i = 1: data.numUEs
    rx(i) = data.ues(i).PhyEntity.StatReceivedPackets;
end
tx = data.gnb.PhyEntity.StatTransmittedPackets;



if(tx>txvalSave)
    data.packet = [data.packet ;[data.Simulator.CurrentTime,zeros(1,data.numUEs)]];
    data.nbsend =data.nbsend + 1;
    txvalSave = tx ;

    %disp("send") 
end

for i = 1:data.numUEs
    if(rx(i)>rxvalSave(i))
        data.packet(ct(i),i+1) = data.Simulator.CurrentTime;
        ct(i) = ct(i)+1 ;
        data.nbreceive(i) = data.nbreceive(i)+ 1;
        rxvalSave(i) = rx(i) ; 
        
        %disp("receive")
    end
end


