#!/usr/bin/env ruby

# This library is used to generate optimization problem
#
# Author:: Matteo Ragni
# Copyright:: Copyright (2015) University of Trento
# License:: Distributes under the same terms as Ruby

gems = %w|nmatrix nmatrix/lapacke colorize pry|
gems.each do |gem|
  #raise LoadError, "Cannot Load Gem :: #{gem}" unless require gem
  require gem
end

##
# Full container for a set of optimization problem completely defined i Ruby.
# Even if it is not convenient in term of computational time. This is for
# sure a very confortable way to prototype algorigthms.
module Optimization
  load './functions.rb'

  ##
  # This class is a general container for the algorithm that I will define
  class Algorithm
    # Objective function
    attr_reader :obj
    # List of constraints for the problem
    attr_reader :cnt
    # Active set for the constraints
    attr_reader :active_set
    # Lagrange function as a proc

    ##
    # Returns the Lagrange function for a defined problem, in the
    def lagrange

    end




    private
    ##
    # This function re-evaluate lagrange procedure each time
    def lagrange=

    end

  end
end
