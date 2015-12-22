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
    attr_reader :objective
    # List of constraints for the problem
    attr_reader :constraints
    # Lagrange function as a proc
    attr_reader :lagrange
    # History of the optimization problem
    attr_reader :history
    # Options for the algorithm
    attr_accessor :options
    # Size for problem
    attr_reader :size
    # Number of equalities and inequalities
    attr_reader :eq_size, :ineq_size

    ##
    # Initialize a new optimization problem
    def initialize(options, objective)
      raise ArgumentError, "options must be an Hash" unless options.is_a? Hash
      raise ArgumentError, "Objective function must be of type ObjectiveFunction" unless
        (objective.is_a? ObjectiveFunction or
         objective.is_a? LinearObjectiveFunction or
         objective.is_a? LinearObjectiveFunctions or
         objective.is_a? QuadraticObjectiveFunction)
      @options = option
      @objective = objective
      @constraints = []
    end

    ##
    # Returns the Lagrange function value for a defined problem, in the form of a `Proc`
    def lagrange(x)
      return @lagrange.call(x)
    end

    private
    ##
    # This function re-evaluate lagrange function (this will keep lagrange call faster)
    def update_lagrange
      @lagrange = Proc.new do |x|
        value = @objective.f(x)
        @contraints.each do |c|
          value += c.lagrange_term(x)
        end
        return value
      end

      @lagrange_x = Proc.new do |x|
        
      end
    end



  end
end
