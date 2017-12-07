classdef PartitionsCombinationIterator
    
    properties
        % a cell array of partitions
        partitions
        indexes
        currentPos = 1
    end
    
    methods
        
        function this = PartitionsCombinationIterator(partitions)
            this.partitions = partitions;
            
            this.indexes = ones(1, length(this.partitions));
        end
        
        function res = getCandidateCount(this, i)
            res = length(this.partitions{i});
        end
        
        function [this, flag, combination] = getNextCombinatrion(this)
            currentPosIndex = this.indexes(this.currentPos);
            if currentPosIndex < this.getCandidateCount(this.currentPos)
                if currentPosIndex < this.getCandidateCount(this.currentPos)
                    this.indexes(this.currentPos) = currentPosIndex + 1;
                else
                    this.currentPos = this.currentPos + 1;

                    % FIXME

                    if this.currentPos > length(this.partitions)
                        flag = false;
                        combination = [];
                        return;
                    end 

                    this.indexes(this.currentPos) = 2;
                end
                
            end
        end
    end
    
