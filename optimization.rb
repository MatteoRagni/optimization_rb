#!/usr/bin/env ruby

# This library is used to generate optimization problem
#
# Author:: Matteo Ragni
# Copyright:: Copyright (2015) University of Trento
# License:: Distributes under the same terms as Ruby

gems = %w|nmatrix nmatrix/lapacke colorize pry-byebug|
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
      @options[:debug] = 0 unless @options[:debug]
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

    ##
    # Returns false if a point is not feasible (does not live in _P_ domain)
    def feasible?(x)
      is = true
      @constraints.each do |c|
        case c.type
        when :equality
          is false if c.f(x) != 0.0
        when :inequality
          is = false if c.f(x) < 0.0
        end
      end
      return is
    end

=begin
    ##
    # Returns the feasible set (all indexes of violated constraints)
    def feasible_set(x)
      ret = []
      @constraints.each do |c|
        next if c.type == :equality

      end
    end
=end

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
    attr_reader :active_set

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
      raise ArgumentError, "Initial guess must be a vector" unless x0.is_a? NMatrix
      raise ArgumentError, "x0 has a wrong size. Must be of size #{@size} x 1" unless (x0.shape == [@size, 1])

      puts "Starting optimization problem".green if @options[:debug] > 0
      puts " - Reordering constraints".green if @options[:debug] > 0
      reorder_constraints
      puts " - Cholesky factorization for H matrix and g vector".green if @options[:debug] > 0
      cholesky

      x = x0

      puts " - Evaluating initial active set".green if @options[:debug] > 0
      active_set(x)
      puts "   A = #{@active_set}".yellow if @options[:debug] > 0

      puts " - Starting solution loop".green if @options[:debug] > 0
      loop do
        puts "  -- Soving quadratic programming problem".green if @options[:debug] > 1
        z, lambdas = q_problem


        if not feasible? z
          x = line_problem(z, x)
          active_set(x)
        else
          x = z
          feasible_set(x, lambdas)
          if @feasible_set == []
            break
          else
            active_set(x)
          end
        end
      end
      return x
    end

    private
    ##
    # Executes a line search to find the best point on the contraints
    def line_problem(z, x)
      t_set = []
      @constraints.each do |c|
        if (c.type == :inequality and c.f(z) < 0.0)
          t = (-c.b - (c.a.transpose.dot(x))[0])/((c.a.transpose.dot((x - z)))[0])
          #raise RuntimeError, "t is too big: 0 <= t = #{t} <= 1" if (t < 0 or t > 1)
          t_set << t
        end
      end
      t_set.each_with_index { |t,i| (@active_set - [i]) + [i] if t == t_set.min }
      ret = x + ((x - z) * t_set.min)
      return ret
    end

    ##
    # Evaluates the active set for the defined step. Also evaluate active inequality set
    def active_set(x)
      @active_set = []
      @constraints.each_with_index do |c, j|
        case c.type
          when :equality
            @active_set << j
          when :inequality
            @active_set << j if c.active? x
        end
      end
    end

    ##
    # Creates a new feasible set
    def feasible_set(x, lambdas)
      @feasible_set = []
      @active_set.each_with_index do |j, i|
        @feasible_set << j if lambdas[i] < 0
      end
    end

    ##
    # Extracts matrix A and vector b given the active set
    def active_set_matrix
      a_ary = []
      b_ary = []
      @active_set.each do |j|
        a_ary << @constraints[j].a.to_flat_array
        b_ary << @constraints[j].b
      end
      a_matrix = NMatrix.new([@size, @active_set.size], a_ary.flatten)
      b_matrix = NMatrix.new([@active_set.size, 1], b_ary.flatten)
      return a_matrix, b_matrix
    end

    ##
    # Evaluate Cholesky decomposition, used in subproblem solution
    def cholesky
      @c, @ct = @objective.s.factorize_cholesky
      @l = @c.invert
      @lt = @l.transpose
      @g = @objective.b
      @hinv = @lt.dot(@lt)
      @d = @hinv.dot(@g)
    end

    ##
    # Move to the first positions equality constraints. Used to define the active set.
    def reorder_constraints
      eq = []
      ineq = []
      @constraints.each do |c|
        eq << c.dup if c.type == :equality
        ineq << c.dup if c.type == :inequality
      end
      @constraints = []
      eq.each do |c|; @constraints << c; end
      ineq.each do |c|; @constraints << c; end
    end

    ##
    # Returns solution for the sub-problem
    def q_problem
      puts "A = " + @active_set.to_s.green
      binding.pry
      if @active_set != []
        a, b = active_set_matrix
        ginv = a.transpose.dot(@hinv.dot(a)).invert

        lambdas = ginv.dot(a.transpose.dot(@d) - b)
        xs      = @hinv.dot(a.dot(lambdas) - @g)
      else
        lambdas = nil
        xs = @d
      end
      puts "x = " + xs.to_s.red
      puts "λ = " + lambdas.to_s.yellow
      return xs, lambdas
    end
  end
end

if $0 == __FILE__ then
  # Definition of the objective function
  m = N[[1.0,2.0,0.0],
        [0.0,3.0,1.0],
        [0.0,1.0,5.0]]
  g = N[[-3.0],
        [-2.0],
        [-1.0]]
  c = 35.0
  obj = Optimization::QuadraticObjectiveFunction.new(m, g, c)
  # The optimum is in x = { 3.66667, -0.66667, 0.33333}

  # Defintion of a constraint on the minimum
  a = NMatrix.new([3,1],[8.0/3.0, -5.0/3.0, -2.0/3.0])
  b = -32.0/3.0
  cnt = Optimization::LinearConstraintFunction.new(:inequality, a, b)

  a = NMatrix.new([3,1],[1.0,0.0,0.0])
  b = 0.0
  cnt_x = Optimization::LinearConstraintFunction.new(:inequality, a, b)

  a = NMatrix.new([3,1],[0.0,-1.0,0.0])
  b = 0.0
  cnt_y = Optimization::LinearConstraintFunction.new(:inequality, a, b)

  a = NMatrix.new([3,1],[0.0,0.0,1.0])
  b = 0.0
  cnt_z = Optimization::LinearConstraintFunction.new(:inequality, a, b)

  # Creation of a new test scenario
  alg = Optimization::QuadraticOptimizer.new({debug: 10 > 0}, obj)
  alg.add_constraint(cnt)
  alg.add_constraint(cnt_x)
  alg.add_constraint(cnt_y)
  alg.add_constraint(cnt_z)
  x = alg.solve(NMatrix.new([3,1], [33.0, 33.0, 33.0]))

  binding.pry
end
