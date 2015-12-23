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
    attr_reader :lagrange, :lagrange_x, :lagrange_xx, :lagrange_lambda
    # History of the optimization problem
    attr_reader :history
    # Options for the algorithm
    attr_accessor :options
    # Size for problem
    attr_reader :size
    # Number of equalities and inequalities
    attr_reader :eq_size, :ineq_size

    ##
    # Initialize a new optimization problem. As for now, multiple dimensions
    # linear functions are not supported
    def initialize(options, objective)
      raise ArgumentError, "options must be an Hash" unless options.is_a? Hash
      raise ArgumentError, "Objective function must be of type ObjectiveFunction" unless
        (objective.is_a? ObjectiveFunction or
         objective.is_a? LinearObjectiveFunction or
         objective.is_a? QuadraticObjectiveFunction)
      @options = option
      @objective = objective
      @size = objective.size
      @constraints = []
      @eq_size = 0
      @ineq_size = 0
    end

    ##
    # Add a new constraint to the problem
    def add_constraint(c)
      raise ArgumentError, "Argument must be a constraint function" unless
        (c.is_a? ConstraintFunction or
         c.is_a? LinearConstraintFunction or
         c.is_a? QuadraticConstraintFunction)
      raise ArgumentError, "Function must accept x of size #{@size}. Is of size #{c.size}" unless @size = c.size
      @eq_size += 1 if c.type == :equality
      @ineq_size += 1 if c.type == :inequality
      @constraints << c
    end

    ##
    # Returns the Lagrange function value for a defined problem, in the form of a `Proc`
    def lagrange(x)
      return @lagrange.call(x)
    end

    ##
    # Returns the Lagrange derivative for a defined problem, in the form of a `Proc`
    def lagrange_x(x)
      return @lagrange_x.call(x)
    end

    ##
    # Returns the Lagrange second derivative for a defined problem, in the form of a `Proc`
    def lagrange_xx(x)
      return @lagrange_xx.call(x)
    end

    ##
    # Returns the Lagrange constraints value for a defined problem, in the form of a `Proc`
    def lagrange_lambda(x)
      return @lagrange_lambda.call(x)
    end

    ##
    # Solves the optimization problem given a starting guess
    def solve(x0)
      raise RuntimeError, "This method is not implemented!"
    end

    private
    ##
    # This function re-evaluate lagrange function (this will keep lagrange call faster)
    def update_lagrange
      @lagrange = Proc.new do |x|
        value = @objective.f(x)
        @contraints.each do |c|
          value -= c.lamda * c.f(x)
        end
        return value
      end

      @lagrange_x = Proc.new do |x|
        value = @objective.g(x)
        @constraints.each do |c|
          value -= c.lambda * c.g(x)
        end
      end

      @lagrange_xx = Proc.new do |x|
        value = @objective.h(x)
        @constraints.each do |c|
          value -= c.lambda * c.h(x)
        end
      end

      @lagrange_lambda = Proc.new do |x|
        value = NMatrix.zeroes([@constraints.size, 1])
        @constraints.each_with_index do |c, i|
          value[i, 0] = c.f(x)
        end
      end
    end
  end
end
