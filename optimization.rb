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
    # Lagrange functions as a proc
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
      @options = options
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
    # Solves the optimization problem given a starting guess
    def solve(x0)
      raise RuntimeError, "This method is not implemented!"
    end

    private
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

  ##
  # Class that represents a quadratic optimizer with active set search.
  class QuadraticOptimizer < Algorithm
    # Active set
    attr_reader :active_set, :active_inequality

    ##
    # Initializer for our quadratic optimizer with active set
    def initialize(options, objective)
      raise ArgumentError, "Objective function must be a QuadraticObjectiveFunction" unless objective.is_a? QuadraticObjectiveFunction
      super(options, objective)
    end

    ##
    # Add a constraint. It could be only of the linear type
    def add_constraint(c)
      raise ArgumentError, "This algorithm accepts only linear constraint" unless c.is_a? LinearConstraintFunction
      super(c)
    end

    ##
    # Solves algorithm when a
    def solve(x0)
      x = x0
      active_set(x)

      @history.x = []
      @history.f = []

      @history.x << x
      @history.f << @objective.f(x)
      loop do
        f, x, lambda = cbqo(x)

        @history.x << x
        @history.f << @objective.f(x)
      end
    end

    private
    ##
    # Evaluates the active set for the defined step. Also evaluate active inequality set
    def active_set(x)
      @active_set = []
      @active_inequality = []
      @constraints.each_with_index do |c, i|
        if c.active? x
          @active_set << i
          @active_inequality << i if c.type == :inequality
        end
      end
    end

    ##
    #
    def cbqo(x)

      return f, x, lambda
    end
  end
end

if $0 == __FILE__ then
  # Definition
  a = N[[1.0,2.0,3.0],
        [4.0,5.0,6.0],
        [7.0,8.0,9.0]]
  b = N[[-30.0],
        [-20.0],
        [-10.0]]
  c = 35.0
  obj = Optimization::QuadraticObjectiveFunction.new(a, b, c)

  alg = Optimization::QuadraticOptimizer.new({}, obj)
end


## Todo
# - aggiungere una funzione di sorting per avere tutte le disuguaglianze in fondo
